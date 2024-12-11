%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Specular Line- Cyst-L11-5v- Plane Wave Transmission
%  Date 22-10-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

%addpath(genpath(pwd));
filepath = pwd;
experiment = 'Cyst_pin_NO_lens';
dataset_owner = 'Hari';
dataset_name = 'Cyst_pin_NO_lens';
Fnum = 1.5;

GPU             = true;   % set to true to simulate using GPU, to false to simulate in CPU
USE_CPU_C_CODE  = false;  % set to true to simulate using CPU cores in the absence of CUDA
SimPlot         = false;   % set to true to view wave propagation.Will be active only if USE_CPU_CORES == false
RecordMovie     = false;  % set to tru to record the movie of simulation else set to false

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
Trans.tx_steering_angles = 51;
Trans.tx_steering_angles_start = 18;
Trans.tx_width = Trans.tx_N_elements * Trans.tx_pitch;% Total width of the transducer [m]
Trans.tx_pw_steering_angles = 0;% [-Trans.tx_steering_angles_start:(2*Trans.tx_steering_angles_start)/(Trans.tx_steering_angles-1):Trans.tx_steering_angles_start];           % Steering angle of the transmit signal [degrees]
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
dx = 0.75e-4;   %0.5e-4;                 % [m]
dy = dx;                      % [m]
dz = dx;                      % [m]

Nx =900 - 2 * PML_X_size;%1600 - 2 * PML_X_size;    % [grid points]
Ny = 700  - 2 * PML_Y_size;%800  - 2 * PML_Y_size;    % [grid points]
%Nz = 16   - 2 * pml_z_size;    % [grid points]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);%, Nz, dz);

% Deifne the axes and probe geometry
x_grid = kgrid.y(:,:);
probe_geometry = linspace(-0.0191, 0.0191, Trans.tx_N_elements)';

% create the time array
t_end = (Nx * dx) * 2.2 / Trans.tx_c0;   % [s]
CFL = 0.09;
kgrid.t_array = makeTime(kgrid, Trans.tx_c0, CFL, t_end);

%% DEFINE THE MEDIUM PARAMETERS
% =========================================================================
% Creating background maps of random scatterers
% =========================================================================
rho0 = 1000;                    % [kg/m^3]
alpha_coeff = 0.05;      % [dB/(MHz^y cm)]
alpha_power = 1.50;
background_map_mean = 1;
background_map_std = 0.01;
background_map = background_map_mean + background_map_std * randn([Nx, Ny]);

% =========================================================================
% Defining medium speed of sound and density
% =========================================================================

medium.sound_speed            = Trans.tx_c0*ones(Nx, Ny).*background_map;
medium.density                = rho0*ones(Nx, Ny).*background_map;
medium.alpha_coeff            = alpha_coeff*ones(Nx, Ny);
medium.alpha_power            = alpha_power;

%% cyst
% cyst-1
x_pos_mat=[25e-3];%Nx*dx/2+10e-3];%[20e-3 : 5e-3 : 40e-3];
y_pos_mat=[Ny*dy/3];%[-10e-3:5e-3:10e-3]%
radius =0.005; 
% cyst-2
x_pos2_mat=[40e-3];%Nx*dx/2+10e-3];%[20e-3 : 5e-3 : 40e-3];
y_pos2_mat=[Ny*dy*2/3];%[-10e-3:5e-3:10e-3]%
radius2 =0.0025;

background_map_mean1 = 1;
background_map_std1 = 0.008;

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx, Ny]);
scattering_c0 = Trans.tx_c0 + 25 + 75 * scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% define properties
sound_speed_map = Trans.tx_c0 * ones(Nx, Ny) .* background_map;
density_map = rho0 * ones(Nx, Ny) .* background_map;

for i=1:length(x_pos_mat)
    scattering_region2 = makeDisc(Nx, Ny, round(x_pos_mat(i,1)/dx), round(y_pos_mat(i,1)/dx), round(radius/dx))...
                       +makeDisc(Nx, Ny, round(x_pos2_mat(i,1)/dx), round(y_pos2_mat(i,1)/dx), round(radius2/dx));
    background_map_S1 = background_map_mean1 + background_map_std1 * randn([size(medium.density(scattering_region2 == 1))]);
    medium.sound_speed(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
    medium.density(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1);
end
%% pins
x_pos_mat=[20e-3:7e-3:50e-3];%Nx*dx/2+10e-3];%[20e-3 : 5e-3 : 40e-3];
y_pos_mat=[Ny*dy/2-15e-3:7e-3:Ny*dy/2+25e-3];%[-10e-3:5e-3:10e-3]%
radius =0.75e-4; 

x_pos_mat = [x_pos_mat.';x_pos_mat(1,3).*ones(5,1)];    % [m]
y_pos_mat = [y_pos_mat(1,3).*ones(5,1);y_pos_mat.'];    % [m]
for i=1:length(x_pos_mat)
   scattering_region2 = makeDisc(Nx, Ny, round(x_pos_mat(i,1)/dx), round(y_pos_mat(i,1)/dx), round(radius/dx));
   medium.sound_speed(scattering_region2 == 1) = 3100; 
   medium.density(scattering_region2 == 1) = 3100/1.5; 
end
%%
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
source_cycles       = 3;        % number of tone burst cycles
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

karray.setArrayPosition([0, 0], 0);
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
            sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        end
        
        sensor_data(:,:,txAng) = karray.combineSensorData(kgrid, sensor_data_grid_points); % Combine the grid points to get in terms of physical transducer elements
        
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
    figure(1),imagesc(probe_geometry.*1000,z_grid(:,1).*1000,20*log10(beamformedDataDASImage));
%       figure(2),imagesc(20*log10(beamformedDataDASImage));
    colormap(gray);
    colorbar;
    vrange = [-60 0];
    caxis(vrange);
    title([dataset_name,'-',num2str(Trans1.tx_Fc(txFreqN)/10^6),'MHz']);
    
    saveas(figure(1),[fsavepath, '\', loadfname, '.fig']);
    
    clearvars dataset;
end

