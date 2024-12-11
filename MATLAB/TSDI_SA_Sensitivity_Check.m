%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check subaperture length influence on TSDI
% Author:Gayathri Malamal
% Year: 2022
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
Fnum = 1.5;
columns = size(x_grid,2);
N_elements = length(probe_geometry);
c0 = 1540;
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

%% Sub-apertures
% Ns_ele = [128,64,48, 32, 16, 1]; % No. of elements for SA sensitivity check
Ns_ele = [64, 32, 16]; % No. of elements for SA sensitivity check

ItApS = 1;
ItApE = 128;
ItApFull = ItApS:ItApE;

%%
% k-Wave p45
% xrange = 105;
% zrange = 396;

% Verasonics Brachioradialis
% xrange = 96;
% zrange = 864;

% [xrange, zrange] = ginput(1);
% xrange = round(xrange);
% zrange = round(zrange);

% if(strcmp(data_type, 'IQ'))
%     pixelS = size(rawDataIQ,1)*(xrange-1) + zrange;
% else
%     pixelS = size(rawDataRF,1)*(xrange-1) + zrange;
% end

region_spec = [];
%roi1 = drawfreehand('Color','r');
roi1 = drawline('SelectedColor','r');%

region_spec = createMask(roi1);
pixelS = find(region_spec(:)>0);
% [zrange, xrange] = find(region_spec > 0);
%%
alphaR = [];
for nrx=1:N_elements
    alphaR(:,nrx) =   -atand( (x_grid(pixelS) - probe_geometry(nrx,1) )./  z_grid(pixelS));
end

focal_delay = zeros(1,length(probe_geometry));

%%
ScatterMatrixPW = [];
% for an1 = 1: length(angle_subset)
%     clc
%     an = an1;
%     if(strcmp(dataSrc, 'kWave'))
%         [ScatterMatrixPW(:,an1),~,~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, w0,Fnum);
%     else
%         [ScatterMatrixPW(:,an1),~,~] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, w0,Fnum);
%     end
% end
% 
% %figure,plot(smooth(abs(ScatterMatrixPW)));
% figure,imagesc(abs(ScatterMatrixPW));

for an1 = 1: length(angle_subset)
    clc
    an = an1;
    if(strcmp(data_type, 'IQ'))
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an1),~,rxApod] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, w0,Fnum);
        else
            [ScatterMatrixPW(:,:,an1),~,rxApod] = scatter_mat_pixel_IQ(acq1,squeeze(rawDataIQ(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, w0,Fnum);
        end
    else
        if(strcmp(dataSrc, 'kWave'))
            [ScatterMatrixPW(:,:,an1),~,rxApod] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,-angle_subset(an),c0, Fnum);
        else
            [ScatterMatrixPW(:,:,an1),~,rxApod] = scatter_mat_pixel_RF(acq1,squeeze(rawDataRF(:,:,aIdx(an))),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,angle_subset(an),c0, Fnum);
        end
    end
end
%%
%figure,plot(smooth(abs(ScatterMatrixPW)));
% if(strcmp(data_type, 'IQ'))
%     ScatterMatrixPW_mean = abs(squeeze(mean(ScatterMatrixPW,1)));
%     ScatterMatrixPW_mean = ScatterMatrixPW_mean./max(ScatterMatrixPW_mean(:));
%     figure,imagesc(angle_subset_deg,1:128,flip(ScatterMatrixPW_mean));
% else
%     ScatterMatrixPW_mean = abs(hilbert(squeeze(mean(ScatterMatrixPW,1))));
%     ScatterMatrixPW_mean = ScatterMatrixPW_mean./max(ScatterMatrixPW_mean(:));
%     figure,imagesc(angle_subset_deg,1:128,ScatterMatrixPW_mean);
% end
% 
% colorbar;
% set(gca,'TickLabelInterpreter','latex')
% xlabel('\bf{Transmit Angle ($ ^o $)}','interpreter','latex');
% ylabel('\bf{Receive Element $\# $ }','interpreter','latex');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',20);
% set(gcf, 'Position',  [100, 100, 400, 400])

if(strcmp(data_type, 'IQ'))
    ScatterMatrixPW_mean = abs(squeeze(mean(ScatterMatrixPW,1)));
    ScatterMatrixPW_mean = ScatterMatrixPW_mean./max(ScatterMatrixPW_mean(:));
else
    ScatterMatrixPW_mean = abs(hilbert(squeeze(mean(ScatterMatrixPW,1))));
    ScatterMatrixPW_mean = ScatterMatrixPW_mean./max(ScatterMatrixPW_mean(:));
end
figure,h=pcolor(angle_subset_deg,1:128,ScatterMatrixPW_mean);
set(h, 'EdgeColor', 'none');
%colorbar;
set(gca,'TickLabelInterpreter','latex')
xlabel('\bf{Transmit Angle ($ ^o $)}','interpreter','latex');
%ylabel('\bf{Receive Angle ($ ^o $) }','interpreter','latex');
ylabel('\bf{Receive Element $\# $ }','interpreter','latex');
set(gca,'YDir','Normal');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',20);
set(gcf, 'Position',  [100, 100, 400, 400])
%%
s = (squeeze(ScatterMatrixPW));
an = 3;
sn = (2*(s(:,an) - min(s(:,an)))./(max(s(:,an)) - min(s(:,an)))) -1;
figure,plot(sn, 'linewidth',2, 'color','k')
set(gca,'TickLabelInterpreter','latex');
xlabel('\bf{Receive Element \#}','interpreter','latex');
ylabel('\bf{Signal Amplitude}','interpreter','latex');
%title('\bf{PW Image}');
set(gca,'fontsize',20);
set(gcf, 'Position',  [100, 100, 400, 400])
axis tight 

%%
alphaROpti = [];

% for jj = 1:length(Ns_ele)
%     envApS = [];
%     alphaR_ApS = [];
%     for i = 1:ItApE-(Ns_ele(jj)-1)
%         ApS = ItApFull(i):ItApFull(i+(Ns_ele(jj)-1));
%         for ii = 1:length(angle_subset)
%             envApS (ii,i) =  abs(sum(ScatterMatrixPW(ApS,ii),1));
%         end
%         center_ApS(i) = ceil(median(ApS));
%         alphaR_ApS(i,:) = squeeze(alphaR(:, center_ApS(i)));
%     end
%     
%     [maxSAVal(jj,:),maxSAIdx(jj,:)] = (max(envApS,[],2));
%     
%     [maxEnergyVal(jj),OptiTxIdx(jj)] = max(max(envApS,[],2),[],1);
%     alphaROpti(jj) = alphaR_ApS(maxSAIdx(jj,OptiTxIdx(jj)));
% end

% for idx = 1:length(alphaROpti)
%     [~, alphaROptiIdx(idx)] = min(abs(alphaROpti(idx)-alphaR));
% end
% 
% alphaROpti
% alphaROptiIdx

uPW = zeros(size(x_grid,1)*size(x_grid,2),length(Ns_ele));
vPW = zeros(size(x_grid,1)*size(x_grid,2),length(Ns_ele));
s1 = zeros(size(x_grid,1)*size(x_grid,2),length(Ns_ele));
c1 = zeros(size(x_grid,1)*size(x_grid,2),length(Ns_ele));

for jj = 1:length(Ns_ele)
    envApS = [];
    alphaR_ApS = [];
    maxSAVal = [];
    maxSAIdx = [];
    maxEnergyVal = [];
    OptiTxIdx = [];
    
    for i = 1:ItApE-(Ns_ele(jj)-1)
        ApS = ItApFull(i):ItApFull(i+(Ns_ele(jj)-1));
        for ii = 1:length(angle_subset)
            if(strcmp(data_type, 'IQ'))
                envApS (:,ii,i) =  abs(sum(ScatterMatrixPW(:,ApS,ii),2));
            else
                envApS (:,ii,i) =  abs(hilbert(sum(ScatterMatrixPW(:,ApS,ii),2)));
            end
        end
        center_ApS(i) = ceil(median(ApS));
        alphaR_ApS(i,:) = squeeze(alphaR(:, center_ApS(i)));
    end
    
    [maxSAVal,maxSAIdx] = (max(envApS,[],3));
    
    [maxEnergyVal,OptiTxIdx] = max(max(envApS,[],3),[],2);
    
    for p = 1:length(pixelS)
        alphaROpti(p,jj) = alphaR_ApS(maxSAIdx(p,OptiTxIdx(p)),p);
    end
    
    s1(pixelS,jj) = sind(alphaROpti(:,jj));
    c1(pixelS,jj) = cosd(alphaROpti(:,jj));
    
    uPW(pixelS,jj) = maxEnergyVal.*s1(pixelS,jj);
    vPW(pixelS,jj) = maxEnergyVal.*c1(pixelS,jj);

end
%%
% kWave 7p6 tilt 45
% z1 = 150:700;
% x1 = 60:128;

% Verasonics needle no tilt IQ
% z1 = 100:320;%1000;
% x1 = 1:80;

% Verasonics Needle no tilt-RF
% z1 = 250:1000;%320:635;
% x1 = 1:60;

% % % Verasonics brachioradialis no visi1 -RF
z1 = 500 :1500;
x1 = 60:128;%69:128;
% 
% % % Verasonics horiz bone 09052022-IQ
% z1 = 150: 350;
% x1 = 1:120;

% % % Verasonics horiz bone 09052022-RF
% z1 = 448:1050;
% x1 = 1:120;

for jj = 1:length(Ns_ele)
    uPW_SA = [];
    vPW_SA = [];
    absMag = [];
    absMagN = [];
    
    uPW_SA = uPW(:,jj);
    vPW_SA = vPW(:,jj);
    
    absMag = sqrt(uPW_SA.^2 + vPW_SA.^2);
    absMagN = absMag./max(max(absMag));
    
%     NonacceptIdx = absMagN>0.2 & absMagN <0.7;
%     uPW_SA(NonacceptIdx)=0;
%     vPW_SA(NonacceptIdx)=0;
    
    uPW_SA = reshape(uPW_SA, size(x_grid,1),size(x_grid,2));
    vPW_SA = reshape(vPW_SA, size(x_grid,1),size(x_grid,2));
    
    figure,
    ax1 = axes;
    imagesc(probe_geometry(x1).*100, z_axis(z1).*100,20*log10(beamformedDataDASPWImage(z1,x1)));
    %imagesc(1:128, z_axis(z1).*100,20*log10(beamformedDataDASPWImage(z1,x1)));
    %figure,imagesc(20*log10(beamformedDataDASPWImage(:,:)));
    
    colormap(gray);
    vrange = [-70 0];
    caxis(vrange);
    shading('interp');
    set(gca,'TickLabelInterpreter','latex')
    xlabel(ax1,'\bf{Width (cm)}','interpreter','latex');
  %  xlabel(ax1,'\bf{Receive Element \#}','interpreter','latex');
    ylabel(ax1,'\bf{Depth (cm)}','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'fontsize',24);
    set(gcf, 'Position',  [100, 100, 600, 600])
    
    zpixel = z1;% 1:size(rawDataPW,1);
    xpixel = x1;% 1:128;
    
    % zpixel = 1:size(z_grid,1);
    % xpixel = 1:128;
    
    %[~, mpidx] = max(absMag(zpixel,xpixel),[],1);
    xp = xpixel(1):1:xpixel(end);
    zp = zpixel(1):1:zpixel(end);%zpixel(mpidx(1:2:end)); %
    %figure
    
    ax2= ax1;
    colormap(ax2,'hot');
    hold on,
         
    quiverC2D(x_grid(zp,xp).*100,z_grid(zp,xp).*100,uPW_SA(zp,xp),-vPW_SA(zp,xp),80,'LineWidth', 4,'Parent',ax2) ;
    set(gca,'YDir','Reverse');
    colormap(ax1,'gray');
    set(gca,'fontsize',24);
    set(gcf, 'Position',  [100, 100, 600, 600])
    
    dateNow = datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM');
%     saveas(gcf,['SA_Sensitivity_',num2str(Ns_ele(jj)),'_',dateNow,'.fig']);
%     saveas(gcf,['SA_Sensitivity_',num2str(Ns_ele(jj)),'_',dateNow,'.svg']);
    %saveas(gcf,['SA_Sensitivity_',num2str(Ns_ele(jj)),'_',dateNow,'.png']);
end

%%
alphaROpti_mean = mean(alphaROpti, 1);
alphaR_mean = mean(alphaR,1);
for idx = 1:length(alphaROpti_mean)
    [~, alphaROptiIdx(idx)] = min(abs(alphaROpti_mean(idx)-alphaR_mean));
end
OptiTxIdx
alphaROpti_mean
alphaROptiIdx
alphaR_mean(alphaROptiIdx)