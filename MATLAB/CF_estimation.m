%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TSDI for a subaperture length of 96: coherence factor estimation
% Author:Gayathri Malamal
% Year: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

selpath1 = 'E:\IITPKD\Research\Research\Jounals_Gayathri\Visualization Tool\Saved Datasets';

[filename, pathname] = uigetfile( ...
    [selpath1,'\*.mat'], ...
    'Pick a dataset');
load(strcat(pathname,filename));

[filename1, pathname1] = uigetfile( ...
    [selpath1,'\*.mat'], ...
    'Pick the beamformed dataset');
load(strcat(pathname1,filename1));

%%
acq1 ='PW';
data_type = 'RF'; %RF/IQ
%data_type = 'IQ'; %RF/IQ

%dataSrc ='kWave';
dataSrc ='Verasonics';

rawDataRF=double(dataset.rawData);

a = dataset.angles;
angle_subset = linspace(min(a),max(a),15);

%% Needle probe tilt
% angle_subset = [linspace(a(1)*180/pi, a(7)*180/pi,7), linspace(a(8)*180/pi, a(73)*180/pi,8) ].*pi/180;

%% MSK tendon
angle_subset = [linspace(a(1)*180/pi, a(68)*180/pi,10), linspace(a(69)*180/pi, a(73)*180/pi,5) ].*pi/180;
%%
angle_subset_deg = angle_subset.*180/pi;

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

z_axis=0.5*1540*timeVector;

[x_grid,z_grid]=meshgrid(probe_geometry,z_axis);
[zeroAng, zeroAngIdx] = min(abs(dataset.angles-0)); %Find zero angle idx
Fnum = 1.5;
columns = size(x_grid,2);
N_elements = length(probe_geometry);

%%
beamformedDataDAS_PW = sum(beamformedDataDASPWTemp(:,:,aIdx), 3);
if(strcmp(data_type, 'IQ'))
    envDAS_PW = abs(beamformedDataDAS_PW);
else
    envDAS_PW = abs(hilbert(beamformedDataDAS_PW));
end

beamformedDataDASPWImage = (envDAS_PW(:,:)./max(max(envDAS_PW(:,:))));
figure,
%ax1 = axes;
%imagesc(probe_geometry.*100, z_axis.*100,20*log10(beamformedDataDASPWImage(:,:)));

imagesc(20*log10(beamformedDataDASPWImage(:,:)));

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
%% Select specular region to test coherence
region_spec = zeros(size(x_grid,1), size(x_grid,2));
roi1 = drawfreehand('Color','r');
%roi1 = drawline('SelectedColor','r');%
region_spec = createMask(roi1);
pixelS = find(region_spec(:)>0);
[zrange, xrange] = find(region_spec > 0);

region_diff = [];
%roi2 = drawfreehand('Color','r');
roi2 = drawline('SelectedColor','b');%
region_diff = createMask(roi2); 
pixelD = find(region_diff(:)>0);

% [xrange, zrange] = ginput(10);
% xrange =  round(xrange);
% zrange =  round(zrange);
% pixelS = (size(rawDataIQ,1).*(xrange(:)-1))+zrange(:);

%% delay compensation

M = 15;%length(angle_subset);
Nc =128;

focal_delay = zeros(1,length(probe_geometry));
c0 = 1540;

%tic;
ScatterMatrixPW = [];
ScatterMatrixPW_Apod =[];
for an1 = 1: length(angle_subset)
    clc
     an = an1; 
     if(strcmp(data_type, 'IQ'))
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an1),ScatterMatrixPW_Apod(:,:,an1),~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, w0,Fnum);
        else
            [ScatterMatrixPW(:,:,an1),ScatterMatrixPW_Apod(:,:,an1),~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, w0,Fnum);
        end
    else
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an1),ScatterMatrixPW_Apod(:,:,an1),~] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, Fnum);
        else
            [ScatterMatrixPW(:,:,an1),ScatterMatrixPW_Apod(:,:,an1),~] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, Fnum);
        end
    end
end

%% Coherence factor calculation
coherent_sum = [];
incoherent_sum =[];

coherent_sum = abs(squeeze(sum(sum(ScatterMatrixPW,2),3))).^2;
incoherent_sum = squeeze(sum(sum(abs(ScatterMatrixPW).^2,2),3));

CF_temp1 = coherent_sum./(M.*Nc.*incoherent_sum);

CF =[ mean(CF_temp1) ]%, mean(CF_temp2)]

coherent_sum_apod = [];
incoherent_sum_apod =[];

coherent_sum_apod = abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3))).^2;
incoherent_sum_apod = squeeze(sum(sum(abs(ScatterMatrixPW_Apod).^2,2),3));

CF_temp2 = coherent_sum_apod./(M.*Nc.*incoherent_sum_apod);

CF_Apod =[ mean(CF_temp2) ]%, mean(CF_temp2)]
%% SNR in dB
mu_spec =  mean(beamformedDataDASPWImage(pixelS)); %mean(abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3))));
std_diff = std(beamformedDataDASPWImage(pixelD)); %std(abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3)))); %
std_spec = std(beamformedDataDASPWImage(pixelS)); %std(abs(squeeze(sum(sum(ScatterMatrixPW_Apod,2),3)))); %

snr_val_db = 20*log10( mu_spec./std_diff)

sig_max = max(beamformedDataDASPWImage(pixelS));
sig_min = min(beamformedDataDASPWImage(pixelS));

%snr_db = 20*log10((sig_max-sig_min)./std_spec)
snr_db = 20*log10(mu_spec./std_spec)

% %% for bone fracture data
% mu_diff =  mean(beamformedDataDASPWImage(pixelD)); 
% snr_db = 20*log10(mu_diff./std_diff)
%% Mean to std ratio (MSD)
meanP = mean(sum(ScatterMatrixPW_Apod,3),2);
stdP = std(sum(ScatterMatrixPW_Apod,3),[],2);
mean_std_R = mean(abs(meanP./stdP))

%% CR in dB
CR_PW = 20*log10(mean(beamformedDataDASPWImage(pixelS))./mean(beamformedDataDASPWImage(pixelD)))
%%
%save (['CFroi_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.mat'], 'CFroi','-v7.3');