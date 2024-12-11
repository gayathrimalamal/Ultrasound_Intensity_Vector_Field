%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TSDI for a subaperture length of 96
%  Author:Gayathri Malamal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

selpath = pwd;

% PW
acq1 = 'PW';
[filename, pathname] = uigetfile( ...
    [selpath,'\*.mat'], ...
    'Pick a dataset');

load(strcat(pathname,filename));

%selpath1 = pwd;
[filename1, pathname1] = uigetfile( ...
    [selpath,'\*.mat'], ...
    'Pick the beamformed dataset');
load(strcat(pathname1,filename1));
%%
data_type = 'RF'; %RF/IQ
%data_type = 'IQ'; %RF/IQ

%dataSrc ='kWave';
dataSrc ='Verasonics';

rawDataRF=double(dataset.rawData);
%rawData=double(dmskp);

%rawData = medfilt1(rawData);
%rawDataRF(:, 106, :) = eps;
a = dataset.angles;
angle_subset = linspace(min(a),max(a),15);

%% Needle probe tilt
%angle_subset = [linspace(a(1)*180/pi, a(7)*180/pi,7), linspace(a(8)*180/pi, a(73)*180/pi,8) ].*pi/180;

%% MSK tendon
%angle_subset = [linspace(a(1)*180/pi, a(68)*180/pi,10), linspace(a(69)*180/pi, a(73)*180/pi,5) ].*pi/180;
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
ItApS = 1;
ItApE = 128;
ItApFull = ItApS:ItApE;
%%
beamformedDataDAS_PW = sum(beamformedDataDASPWTemp(:,:,:), 3);
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

%% kWave 7p6 tilt 0
z1 = 200:600;
x1 = 35:100;
% kWave 7p6 tilt 20
% z1 = 200:700;
% x1 = 50:128;
% % kWave 7p6 tilt 45
% z1 = 150:700;
% x1 = 60:128;
% 
% % kWave 7p6 probe tilt 45
% z1 = 100:500;
% x1 = 65:120;

%% Verasonics: Define ROI 

% bone fracture1
% z1 = 75: 250;
% x1 = 10:105;
% % Verasonics horizontal bone
% z1 = 150: 350;
% x1 = 30:90;
% Verasonics needle no tilt-IQ
z1 = 100: 320;
x1 = 1:60;

% Verasonics Needle no tilt-RF
% z1 = 250:1000;%320:635;
% x1 = 1:60;

%Verasonics Needle tube tilt
% z1 = 250:1000;%320:635;
% x1 = 10:100;%20:120;%

% % Verasonics horiz bone 09052022
% z1 = 150: 350;
% x1 = 1:128;
% % Verasonics cross bone 09052022
% z1 = 50: 200;
% x1 = 15:90;
% % Verasonics brachioradialis no visi1 -RF
% z1 = 500 :1500;
% x1 = 60:128;%69:128;

% Verasonics brachioradialis no visi1 -IQ
% z1 = 168:500;
% x1 = 60:128;%69:128;

% Verasonics brachioradialis no visi1 new -RF
z1 = 100:1000;
x1 = 40:128;%69:128;
%%
% figure,
% ax1 = axes;
% imagesc(probe_geometry(x1).*100, z_axis(z1).*100,20*log10(beamformedDataDASPWImage(z1,x1)));
% 
% %figure,imagesc(20*log10(beamformedDataDASPWImage(:,:)));
% 
% colormap(gray);
% %colorbar;
% vrange = [-70 0];
% caxis(vrange);
% shading('interp');
% set(gca,'TickLabelInterpreter','latex')
% colorbar('TickLabelInterpreter','latex');
% hold on;
% xlabel(ax1,'\bf{x (channels)}','interpreter','latex');
% ylabel(ax1,'\bf{z (time samples)}','interpreter','latex');
% %title('\bf{PW Image}');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',24);
% set(gcf, 'Position',  [100, 100, 600, 600])
% hold off;
%%
tic
% Select the ROI on the original figure and not on the reduced size figures as it will change the sample indices
region_spec = zeros(size(x_grid,1), size(x_grid,2));
roi1 = drawfreehand('Color','r');
%roi1 = drawline('SelectedColor','r');%
region_spec = createMask(roi1);
pixelS = find(region_spec(:)>0);
[zrange, xrange] = find(region_spec > 0);

%% If selecting two ROIs

% Select the ROI on the original figure and not on the reduced size figures as it will change the sample indices
% region_spec1 = zeros(size(x_grid,1), size(x_grid,2));
% %roi1 = drawfreehand('Color','r');
% roi1 = drawline('SelectedColor','r');%
% region_spec1 = createMask(roi1);
% 
% roi2 = drawline('SelectedColor','r');%
% region_spec2 = createMask(roi2);
% 
% region_spec = region_spec1+ region_spec2;
% pixelS = find(region_spec(:)>0);
% [zrange, xrange] = find(region_spec > 0);

%%
alphaR = [];
for nrx=1:N_elements
    alphaR(:,nrx) =   -atand( (x_grid(pixelS) - probe_geometry(nrx,1) )./  z_grid(pixelS));
end

s1 = zeros(size(x_grid,1)*size(x_grid,2),1);
c1 = zeros(size(x_grid,1)*size(x_grid,2),1);

focal_delay = zeros(1,length(probe_geometry));
c0 = 1540;

%tic;
%parpool;
ScatterMatrixPW = [];
for an = 1: length(angle_subset)
    clc
    an
    if(strcmp(data_type, 'IQ'))
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an),~,~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, w0,Fnum);
        else
            [ScatterMatrixPW(:,:,an),~,~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, w0,Fnum);
        end
    else
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an),~,~] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, Fnum);
        else
            [ScatterMatrixPW(:,:,an),~,~] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, Fnum);
        end
    end
end

envApS = [];
alphaR_ApS = [];
alphaROpti = [];
stdP =[];
meanP = [];
mstd_R = zeros(size(x_grid,1)*size(x_grid,2),1);

for i = 1:ItApE-31
    ApS = ItApFull(i):ItApFull(i+31);
    for ii = 1:length(angle_subset)
        envApS (:,ii,i) =  abs(sum(ScatterMatrixPW(:,ApS,ii),2));
    end
    center_ApS(i) = ceil(median(ApS));
    alphaR_ApS(i,:) = squeeze(alphaR(:, center_ApS(i)));
end

[maxSAVal,maxSAIdx] = (max(envApS,[],3));

[maxEnergyVal,OptiTxIdx] = max(max(envApS,[],3),[],2);

for p = 1:length(pixelS)
    alphaROpti(p) = alphaR_ApS(maxSAIdx(p,OptiTxIdx(p)),p);
    meanP(p) = mean(ScatterMatrixPW(p,:,OptiTxIdx(p)));
    stdP(p) = std(ScatterMatrixPW(p,:,OptiTxIdx(p)));
end
mstd_R(pixelS) = abs(meanP./stdP);

s1(pixelS) = sind(alphaROpti(:));
c1(pixelS) = cosd(alphaROpti(:));
%%
uPW = zeros(size(x_grid,1)*size(x_grid,2),1);
vPW = zeros(size(x_grid,1)*size(x_grid,2),1);
uPW(pixelS) = maxEnergyVal.*s1(pixelS);
vPW(pixelS) = maxEnergyVal.*c1(pixelS);

% uPW(find(mstd_R>0.35)) = 0;
% vPW(find(mstd_R>0.35)) = 0;

uPW = reshape(uPW, size(x_grid,1),size(x_grid,2));
vPW = reshape(vPW, size(x_grid,1),size(x_grid,2));

%%
figure,
ax1 = axes;
imagesc(probe_geometry(x1).*100, z_axis(z1).*100,20*log10(beamformedDataDASPWImage(z1,x1)));

%figure,imagesc(20*log10(beamformedDataDASPWImage(:,:)));

colormap(gray);
vrange = [-70 0];
caxis(vrange);
shading('interp');
set(gca,'TickLabelInterpreter','latex')
%colorbar('TickLabelInterpreter','latex');
xlabel(ax1,'\bf{Width (cm)}','interpreter','latex');
ylabel(ax1,'\bf{Depth (cm)}','interpreter','latex');
%title('\bf{PW Image}');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',24);
set(gcf, 'Position',  [100, 100, 600, 600])
%%
% saveas(gcf,['PW_',num2str(length(angle_subset)),'_DAS','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.fig']);
% saveas(gcf,['PW_',num2str(length(angle_subset)),'_DAS','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.svg']);
%%
zpixel = z1;% 1:size(rawDataPW,1); 
xpixel = x1;% 1:128; 

%figure
xp = xpixel(1):3:xpixel(end);%xpixel(1):2:xpixel(end);
zp = zpixel(1):7:zpixel(end);%zpixel(mpidx(1:2:end)); %
ax2= ax1;
colormap(ax2,'hot');
hold on,

uPWs = uPW; %smoothn(uPW, 'robust');%
vPWs = vPW; %smoothn(vPW, 'robust');%

% uPWs = smoothn(uPW, 'robust');%
% vPWs = smoothn(vPW, 'robust');%
%quiver(x_grid(zp,xp).*100,z_grid(zp,xp).*100,uPWs(zp,xp),-vPWs(zp,xp),7) ;
quiverC2D(x_grid(zp,xp).*100,z_grid(zp,xp).*100,uPWs(zp,xp),-vPWs(zp,xp),7,'LineWidth', 3,'Parent',ax2) ;

%axis tight
set(gca,'YDir','Reverse');
%axis equal
%linkaxes([ax1,ax2]);
%ax2.Visible = 'off';
colormap(ax1,'gray');
%colormap(ax2,'hot');
set(gca,'fontsize',24);
set(gcf, 'Position',  [100, 100, 600, 600])
%colorbar(ax2);%'Position',[.88 .11 .0175 .815],'TickLabelInterpreter','latex','FontWeight','Bold');
%colorbar(ax2,'Position',[.92 .11 .0175 .815],'TickLabelInterpreter','latex','FontWeight','Bold');
toc
%%
saveas(gcf,['PW',num2str(length(angle_subset)),'DASQuiverCOverlayRoi','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.fig']);
saveas(gcf,['PW',num2str(length(angle_subset)),'DASQuiverCOverlayRoi','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.svg']);
