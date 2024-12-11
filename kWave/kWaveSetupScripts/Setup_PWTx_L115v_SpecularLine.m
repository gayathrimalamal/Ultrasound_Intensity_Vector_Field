%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Specular Line- Specular Reflector Model-L11-5v- Plane Wave Transmission
%  Date 22-10-2021
%  Author: Gayathri Malamal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%close all;

%addpath(genpath(pwd));
filepath = pwd;
experiment = 'specularLine';
dataset_owner = 'Gayathri';
dataset_name = 'specularLine';

GPU             = false;   % set to true to simulate using GPU, to false to simulate in CPU
USE_CPU_C_CODE  = false;  % set to true to simulate using CPU cores in the absence of CUDA
SimPlot         = true;   % set to true to view wave propagation.Will be active only if USE_CPU_CORES == false
RecordMovie     = true;  % set to tru to record the movie of simulation else set to false

%% Define L11-5V transducer
%endDepth = 3e-2;                      % Maximun depth in diagonal dir [m]
Trans.name = 'L11-5v';                % Transducer name
Trans.tx_N_elements = 128;            % Number of transducer elements
Trans.tx_element_width  = 270e-6;      % Width of a single transducer element [m]
Trans.tx_kerf = 30e-6;                % Kerf of the transducer [m]
Trans.tx_pitch = 300e-6;              % Pitch of the transducer [m]
Trans.tx_elevation_height = 5e-3;     % Elevation height of the transducer [m]
Trans.tx_elevation_focus = 18e-3;     % Elevation focus of the transducer [m]
Trans.tx_Fc = 7.6e6;   %2e6;          % Transmit frequencies [Hz]
Trans.tx_Fmax = 15e6;                 % Maximum possible transducer center frequency [Hz]
Trans.tx_c0 = 1540;                   % SoS of the transducer [m/s]
%Trans.tx_steering_angles = 51;
Trans.tx_steering_angles_start = 18;
Trans.tx_width = Trans.tx_N_elements * Trans.tx_pitch;% Total width of the transducer [m]
Trans.tx_pw_steering_angles = 0;%linspace(-18,18,73);% [-Trans.tx_steering_angles_start:(2*Trans.tx_steering_angles_start)/(Trans.tx_steering_angles-1):Trans.tx_steering_angles_start];           % Steering angle of the transmit signal [degrees]
Trans.na = length(Trans.tx_pw_steering_angles);       % Number of steering angles
Trans.tx_active_SA = 1: Trans.tx_N_elements;          % Number of active transducer elements at transmit
Trans.rx_active_SA = 128;                             % Number of active transducer elements at a receive
rx_SubAperture = 128;                                 % Subaperture widths to be extracted separately from rx_active

%% DEFINE THE K-WAVE GRID

%Important note: In k-wave x is the depth direction, y is the sensor (or lateral) direction and z is the elevation

% set the size of the perfectly matched layer (PML)
PML_X_size = 10;
PML_Y_size = 10; 

% set desired grid size in the x-direction not including the PML
% x = 0.06;                      % [m] % This is the depth of the lung-->6cm

c_min = Trans.tx_c0;
% lambda = c_min/Trans.tx_Fc ;
lambda = c_min/(7.6e6) ;

% Calculate the spacing between the grid points
dx = 0.75e-4;  %lambda/2;       % [m]
dy = dx;                      % [m]
dz = dx;                      % [m]

Nx = 300 - 2 * PML_X_size;%1600 - 2 * PML_X_size;    % [grid points]
Ny = 534  - 2 * PML_Y_size;%800  - 2 * PML_Y_size;    % [grid points]
%Nz = 16   - 2 * pml_z_size;    % [grid points]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);%, Nz, dz);

% Deifne the axes and probe geometry
x_grid = kgrid.y(:,:);
probe_geometry = linspace(-0.0191, 0.0191, Trans.tx_N_elements)';

% create the time array
t_end = (Nx * dx) * 2.2 / Trans.tx_c0;   % [s]
CFL = 0.05;
kgrid.t_array = makeTime(kgrid, Trans.tx_c0, CFL, t_end);

%% DEFINE THE MEDIUM PARAMETERS
% =========================================================================
% Creating background maps of random scatterers
% =========================================================================
rng(1);
scat_dist = randn(Nx,Ny);
background_map_mean = 1;
background_map_std = 0.01;
background_map = background_map_mean + background_map_std * randn([Nx, Ny]);

% =========================================================================
% Defining medium speed of sound and density
% =========================================================================
% % -20 tilt
% lineStartX = round(Nx/2)+10; %Start X point of the line
% lineStartY = round(Ny/2)-50; %Start Y point of the line

% 0 tilt
lineStartX = round(Nx/2)-30; %Start X point of the line
lineStartY = round(3*Ny/4)-20;%round(Ny/2)-50; %Start Y point of the line

% % +20 tilt
% lineStartX = round(Nx/2)-60; %Start X point of the line
% lineStartY = round(3*Ny/4)+80; %Start Y point of the line
% 
% % +45 tilt
lineStartX = round(Nx/2)-90; %Start X point of the line
lineStartY = round(3*Ny/4)+90; %Start Y point of the line

lineLen = 200;%Length of the line
lineThickness = 10; %Thickness of the line
tilt =45 ; %Tilt of the line

[~, sound_speed, density, alpha_coeff, Fnum] = specularReflLine1(Nx, Ny, lineStartX, lineStartY, lineLen, lineThickness, tilt, scat_dist);
medium.sound_speed            = sound_speed;
medium.density                = density;
medium.alpha_coeff            = alpha_coeff;
medium.alpha_power            = 1.5;
medium.BonA = 6;
%medium.alpha_mode = 'no_dispersion';

figure,imagesc(medium.sound_speed); % Plotting the speed of sound map
title('SoS map');
colorbar;
figure,imagesc(medium.density); % Plotting the density map
title('Density map');
colorbar;
%% DEFINE THE INPUT SIGNAL
% define properties of the input signal
% Define the driving signal
source_strength     = 4e6;      % [Pa], Maximum is 4e6
source_cycles       = 4;        % number of tone burst cycles
source.p_mask   = zeros(Nx, Ny);
source_width    = round(Trans.tx_element_width/dy);
source_kerf     = round(Trans.tx_kerf/dy);
element_height  = round(Trans.tx_elevation_height/dx);
rotation        = 0;

% Create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% Creating each transducer element by combining grid points
for Nc = 1: Trans.tx_N_elements
    karray.addRectElement([kgrid.x_vec(1), probe_geometry(Nc, 1)],100e-6, Trans.tx_element_width, rotation);
end

karray.setArrayPosition([0, 0], rotation);
source.p_mask = karray.getArrayBinaryMask(kgrid);
%source.p_mode='dirichlet';

display_mask=medium.sound_speed/norm(medium.sound_speed);
display_mask=histeq(display_mask);
display_mask(display_mask<0.7)=0;
display_mask(display_mask>0.7)=1;

sensor.mask = zeros(Nx, Ny);
sensor.mask = source.p_mask;

%% RUN THE SIMULATION

for txFreqN = 1: length(Trans.tx_Fc)
    acq = 'PW';
    source_freq         = Trans.tx_Fc(txFreqN);    % [Hz]
    for txAng = 1: length(Trans.tx_pw_steering_angles)
        
        disp(txAng);
        
        % Calculate the steering delay
        if(Trans.tx_pw_steering_angles(txAng)>0)
            tx_steering_delay(:,txAng) = (probe_geometry(Trans.tx_N_elements,1)-probe_geometry(Trans.tx_active_SA,1))*sind(Trans.tx_pw_steering_angles(txAng))/(Trans.tx_c0*kgrid.dt);
        else
            tx_steering_delay(:,txAng) = (probe_geometry(1,1)-probe_geometry(Trans.tx_active_SA,1))*sind(((Trans.tx_pw_steering_angles(txAng))))/(Trans.tx_c0*kgrid.dt);
        end
        
        % Note that the source strength is not scaled by acoustic impedance as a pressure source is used (Z=P/V)
        source_sig = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles,'SignalOffset',tx_steering_delay(:,txAng));
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        
        if(GPU)
            DATA_CAST       = 'gpuArray-single';
            input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PlotPML', false, ...
                'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
                'DisplayMask', display_mask, 'DataCast', DATA_CAST, 'DataRecast', true,...
                'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'}; % Video of the simulation will be recorded in the current directory if set to true
            if(USE_CPU_C_CODE == true)
                sensor_data_grid_points = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
            elseif(SimPlot == true)
                sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
            else
                sensor_data_grid_points = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
            end
        else
            DATA_CAST       = 'single';
            if(SimPlot == true)
                input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PlotPML', false, ...
                    'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
                    'DisplayMask', display_mask, 'DataCast', DATA_CAST, ...
                    'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'}; % Video of the simulation will be recorded in the current directory if set to true
            else
                input_args = {'PMLSize', [PML_X_size, PML_Y_size], 'PMLInside', false, 'DataCast', DATA_CAST,  'DataRecast', true,'PlotSim', false};
            end
            sensor_data_grid_points= kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        end
        sensor_data(:,:,txAng) = karray.combineSensorData(kgrid, sensor_data_grid_points); % Combine the grid points to get in terms of physical transducer elements
%        sensor_data_grid_points1 = karray.combineSensorData(kgrid, sensor_data_grid_points); % Combine the grid points to get in terms of physical transducer elements
%         for idx = 1:128
%             idx
%             sensor_data(:,idx,txAng) = attenComp(sensor_data_grid_points1(idx,:), kgrid.dt, medium.sound_speed, 0.5,1.5);
%         end
        
        
    end
    RcvData = permute(sensor_data, [2 1 3]);
    dataset = dataExtract_MultiFc_kWave(filepath, experiment, dataset_name, dataset_owner, acq, Trans, tx_steering_delay, probe_geometry, RcvData, kgrid, medium,txFreqN);
end

%
%% DAS Beamforming

%clearvars dataset

% DAS PW beamforming
acq = 'PW';
if ~exist([filepath,'\kWaveBeamformedImages\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name], 'dir')
    mkdir([filepath,'\kWaveBeamformedImages\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name]);
end

fsavepath = [filepath,'\kWaveBeamformedImages\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name];

for txFreqN = 1: length(Trans.tx_Fc)
    txFreqN
    [loadfilename, loadpath] = uigetfile( ...
        [filepath,'\kWaveDatasets\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name,'\*',num2str(Trans.na),'*',dataset_name,'*',num2str(Trans.tx_Fc(txFreqN)/10^6),'*.mat'], ...
        'Pick a dataset to be beamformed');
    
    d=load(strcat(loadpath, loadfilename));
    [loadpath1,loadfname,loadfileext] = fileparts(loadfilename);
    
%     timeVector = kgrid.t_array;
    rawData=double(d.dataset.rawData);
    Fs = d.dataset.Fs*1e6;
    timeVector = (0:(size(rawData,1)-1)).'/(Fs);
    
    z_axis = 0.5*1540*timeVector;
    Trans1 = d.dataset.Trans;
    probe_geometry1 = d.dataset.probe_geometry;
    
    %[x_grid,z_grid] = meshgrid(probe_geometry,z_axis);
    z_axis1 = z_axis;%linspace(min(z_axis),max(z_axis), Nx);
    [x_grid,z_grid] = meshgrid(probe_geometry,z_axis1);
    
    ang = d.dataset.angles;
      
    for txAng = 1: length(ang)
        beamformedDataPW (:,:,txAng) = DAS_column_PW(squeeze(rawData(:, :, txAng)), timeVector, x_grid, z_grid, probe_geometry1, -ang(txAng), Trans1.tx_c0*ones(size(x_grid)), Fnum);
    end
    
    beamformedDataDAS = sum(beamformedDataPW, 3);
    env = abs(hilbert(beamformedDataDAS));
    beamformedDataDASImage = (env(:,:)./max(max(env(:,:))));
    figure,imagesc(probe_geometry.*1000,z_grid(:,1).*1000,20*log10(beamformedDataDASImage));
    colormap(gray);
    colorbar;
    vrange = [-60 0];
    caxis(vrange);
    title([dataset_name,'-',num2str(Trans1.tx_Fc(txFreqN)/10^6),'MHz']);
    
    saveas(figure(1),[fsavepath, '\', loadfname, '.fig']);
    
    clearvars dataset;
end

