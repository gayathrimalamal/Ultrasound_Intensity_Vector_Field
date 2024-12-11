%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TSDI for a subaperture length of 96: contrast ratio estimation
% Author:Gayathri Malamal
% Year: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clearvars;

%load dataset
selpath = pwd;
dataset_owner = 'Gayathri';
filepath1 = pwd;

[loadfilename, loadpath] = uigetfile( ...
    [filepath1, dataset_owner,'\*','.mat'], ...
    'Pick PW RF dataset');
load(strcat(loadpath, loadfilename));

filepath3 = 'E:\IITPKD\Research\Research\Jounals_Gayathri\Study of Transmit Schemes\';
[loadfilename1, loadpath1] = uigetfile( ...
    [filepath3, dataset_owner,'\*','.mat'], ...
    'Pick PW IQ bf data');
load(strcat(loadpath1, loadfilename1));
%%
acq2 = 'PW';
dataSrc = 'kWave';
%dataSrc = 'Verasonics';
%data_type = 'RF'; %RF/IQ
data_type = 'IQ'; %RF/IQ

rawDataPW = double(dataset.rawData);

angles = dataset.angles;
angle_subset = linspace(min(a),max(a),15);
%% Needle probe tilt
%angle_subset = [linspace(a(1)*180/pi, a(7)*180/pi,7), linspace(a(8)*180/pi, a(73)*180/pi,8) ].*pi/180;

%% MSK tendon
%angle_subset = [linspace(a(1)*180/pi, a(68)*180/pi,10), linspace(a(69)*180/pi, a(73)*180/pi,5) ].*pi/180;
%%
if(strcmp(data_type, 'IQ'))
    
    if(strcmp(dataSrc, 'Verasonics'))
        [~,aIdx] = min(abs(a-angle_subset));
        probe_geometry = dataset.probe_geometry(:,1)*1540/(dataset.Trans.frequency*1e6);
        for i = 1: size(rawDataRF, 3)
            [rawDataIQ(:,:,i),Fc,Fs] = rf2iqFn(rawDataRF(:,:,i), dataset.Fs*1e6, dataset.Trans.frequency*1e6);
        end
    else
        probe_geometry = dataset.probe_geometry(:,1);
        [~,aIdx] = min(abs(a'-angle_subset));
        for i = 1:  size(rawDataRF, 3)
            [rawDataIQ(:,:,i),Fc,Fs] = rf2iqFnkWave(rawDataRF(:,:,i), dataset.Fs*1e6, dataset.Fc*1e6);
        end
    end
    w0 = 2*pi*Fc;
    timeVector=(0:size(rawDataIQ,1)-1)'/(Fs);
else
    if(strcmp(dataSrc, 'Verasonics'))
        [~,aIdx] = min(abs(a-angle_subset));
        probe_geometry = dataset.probe_geometry(:,1)*1540/(dataset.Trans.frequency*1e6);
    else
         probe_geometry = dataset.probe_geometry(:,1);
        [~,aIdx] = min(abs(a'-angle_subset));
    end
    timeVector=(0:size(rawDataRF,1)-1)'/(dataset.Fs*1e6);
end

%%
c0 = dataset.c0;
z_axis = 0.5*1540*timeVector;
[x_grid,z_grid] = meshgrid(probe_geometry,z_axis);
focal_delay = zeros(1,length(probe_geometry));
Fnum = 1.5;

%%
beamformedDataDAS_PW = sum(beamformedDataDASPWTemp(:,:,aIdx), 3);
envDAS_PW = abs(beamformedDataDAS_PW);
beamformedDataDASPWImage = (envDAS_PW(:,:)./max(max(envDAS_PW(:,:))));
figure,imagesc(20*log10(beamformedDataDASPWImage(:,:)));

colormap(gray);
colorbar;
vrange = [-70 0];
caxis(vrange);
shading('interp');
set(gca,'TickLabelInterpreter','latex')
colorbar('TickLabelInterpreter','latex');
hold on;
xlabel('\bf{x (channels)}','interpreter','latex');
ylabel('\bf{z (time samples)}','interpreter','latex');
title('\bf{PW Image}');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',24);
set(gcf, 'Position',  [100, 100, 600, 600])
hold off;
%%
region_spec = [];
%roi1 = drawfreehand('Color','r');
roi1 = drawline('SelectedColor','r');%
region_spec = createMask(roi1); 
pixelS = find(region_spec(:)>0);

region_diff = [];
%roi2 = drawfreehand('Color','r');
roi2 = drawline('SelectedColor','b');%
region_diff = createMask(roi2); 
pixelD = find(region_diff(:)>0);
%% Contrast Ratio dB
CR_PW = 20*log10(mean(beamformedDataDASPWImage(pixelS))./mean(beamformedDataDASPWImage(pixelD)))
%% SNR in dB
mu =  mean(beamformedDataDASPWImage(pixelS)); %mean(abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3))));
sig = std(beamformedDataDASPWImage(pixelD)); %std(abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3)))); %
snr_val_db = 20*log10( mu./sig)