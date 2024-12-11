%% To extract necessary parameters and variables for beamforming
function [envr]=Reconstruct_Vector_Field(rdata)
Trans = evalin('base','Trans');
TX = evalin('base','TX');
TW = evalin('base','TW');
Resource = evalin('base','Resource');
P = evalin('base','P');
PData = evalin('base','PData');
na = evalin('base','na');
dtheta = evalin('base','dtheta');
Receive = evalin('base','Receive');
% runAcqPeriod = evalin('base','ImgData');

startDepth = P.startDepth* Resource.Parameters.speedOfSound/(Trans.frequency*1e6);
endDepth = P.endDepth* Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

dataRF= rdata(Receive(1).startSample : Receive(na).endSample, :,1);
dataRF = reshape(dataRF,[Receive(1).endSample, na, size(dataRF, 2)]);
dataRF = permute(dataRF, [1, 3, 2]);

angI = round(linspace(1,73,15));

for angleIndex = 1:length(angI)
    delay = TX(angI(angleIndex)).Delay/(Trans.frequency*1e6);
    maxDelay = (max(delay)/2)*(Receive(1).decimSampleRate*1e6);
    
    if(maxDelay==0)
        maxDelay=1;
    end
    
    rawData1 = squeeze(dataRF(round(maxDelay):end,:,angI(angleIndex)));
    rawDataRF(1:size(rawData1,1),:,angleIndex) = double(rawData1);
end

for i = 1:length(angI)
 %   [rawDataIQ1(:,:,i),~,Fs] = rf2iqFn(rawDataRF(:,:,i), Receive(1).decimSampleRate*1e6, Trans.frequency*1e6);
     [rawData(:,:,i),~,Fs] = rf2iqFn(rawDataRF(:,:,i), Receive(1).decimSampleRate*1e6, Trans.frequency*1e6);
end

%save(['dataset.mat'], 'dataset','-v7.3');
%%
acq = 'PW';
%rawData=double(rawDataIQ1);
timeVector=(0:size(rawData,1)-1)'/(Fs);

if(strcmp(Trans.units, 'mm'))
    probe_geometry = Trans.ElementPos(:,1);
else
    probe_geometry = Trans.ElementPos(:,1)*Resource.Parameters.speedOfSound/(Trans.frequency*1e6);
end

z_axis_IQ = 0.5*1540*timeVector(1:3:end);
z_axis_IQ = z_axis_IQ(z_axis_IQ <0.025);
[x_grid,z_grid]=meshgrid(probe_geometry,z_axis_IQ);
Fnum = 1.5;
w0 = 2*pi*(Receive(1).demodFrequency)*1e6;
c0 = 1540;
N_elements = length(probe_geometry);

%%
focal_delay = zeros(1,size(rawData,2));
for txAng=1:length(angI)
    %txAng
  AngIdx = txAng;
  Angl = TX(angI(txAng)).Steer(1);%P.startAngle+(txAng-1)*dtheta);
  beamformedDataDASPW(:,:,txAng) = DAS_pixel_IQ(acq,rawData(:,:,AngIdx),timeVector,focal_delay, x_grid,z_grid,probe_geometry,Angl,c0*ones(size(x_grid)), w0,Fnum);
end

z_axis_r = linspace(startDepth,endDepth,PData.Size(1));

beamformedDataDAS = sum(beamformedDataDASPW, 3);
envn =  abs(beamformedDataDAS);
envn = envn./max(envn(:));
env = abs(beamformedDataDAS).*10^3;
envr = interp1(z_axis_IQ,env,z_axis_r,'linear',0);

figure(8);
I = 20*log10(envn);
imagesc(probe_geometry.*1000, z_axis_IQ.*1000,I);

colormap(gray);
%colorbar;
vrange = [-70 0];
caxis(vrange);
shading('interp');
%set(gca,'TickLabelInterpreter','latex')
%colorbar('TickLabelInterpreter','latex');
hold on;
xlabel('\bf{mm}')
% xlabel('\bf{x (channels)}','interpreter','latex');
% ylabel('\bf{z (time samples)}','interpreter','latex');
% title('\bf{PW Image}');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',24);
%set(gcf, 'Position',  [100, 100, 600, 600])
%hold off;
%%
figure(9);
ax1 = axes;
I = 20*log10(envn);
imagesc(probe_geometry.*1000, z_axis_IQ.*1000,I);

colormap(gray);
%colorbar;
vrange = [-70 0];
caxis(vrange);
shading('interp');
%set(gca,'TickLabelInterpreter','latex')
%colorbar('TickLabelInterpreter','latex');
hold on;
xlabel(ax1,'\bf{mm}')
% xlabel('\bf{x (channels)}','interpreter','latex');
% ylabel('\bf{z (time samples)}','interpreter','latex');
% title('\bf{PW Image}');
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',24);
%set(gcf, 'Position',  [100, 100, 600, 600])
%hold off;

%% TSDI Estimation
ItApS = 1;
ItApE = 128;
ItApFull = ItApS:ItApE;
threshold = 0.15; % between 0 and 1
% Left ROI
% % Bounding Box Dimensions in mm
minX = probe_geometry(5).*1000;
maxX = probe_geometry(64).*1000;

% Center ROI
% Bounding Box Dimensions in mm
% minX = probe_geometry(32).*1000;
% maxX = probe_geometry(96).*1000;

% % Right ROI
% % Bounding Box Dimensions in mm
% minX = probe_geometry(64).*1000;
% maxX = probe_geometry(126).*1000;
% 
% %Full ROI 
% % Bounding Box Dimensions in mm
% minX = probe_geometry(5).*1000;
% maxX = probe_geometry(126).*1000;

minZ = 13;
maxZ = 23;

region_spec = (z_grid>(minZ*1e-3)) & (z_grid<(maxZ*1e-3)) & (x_grid>(minX*1e-3)) & (x_grid<(maxX*1e-3));
region_spec(1:2:end,:) = zeros(size(region_spec(1:2:end,:)));
region_spec(:,1:2:end) = zeros(size(region_spec(:,1:2:end)));
pixelS = find(region_spec(:)>0);

region_full = (z_grid>(minZ*1e-3)) & (z_grid<(maxZ*1e-3)) & (x_grid>(minX*1e-3)) & (x_grid<(maxX*1e-3));%zeros(size(x_grid,1), size(x_grid,2));
%region_full(:,1:128) = 1;
pixelFull = find(region_full(:)>0);

figure(9),

plot([minX maxX],[minZ minZ],'g-','linewidth',2);
plot([minX maxX],[maxZ maxZ],'g-','linewidth',2);
plot([minX minX],[minZ maxZ],'g-','linewidth',2);
plot([maxX maxX],[minZ maxZ],'g-','linewidth',2);
hold off;

%%
alphaR = [];
for nrx=1:N_elements
    alphaR(:,nrx) =   -atand( (x_grid(pixelS) - probe_geometry(nrx,1) )./  z_grid(pixelS));
end

uPW = zeros(size(x_grid,1)*size(x_grid,2),1);
vPW = zeros(size(x_grid,1)*size(x_grid,2),1);
s1 = zeros(size(x_grid,1)*size(x_grid,2),1);
c1 = zeros(size(x_grid,1)*size(x_grid,2),1);

ScatterMatrixPW = [];
for an = 1: length(angI)
     Angl = TX(angI(txAng)).Steer(1);
    [ScatterMatrixPW(:,:,an),~,~] = scatter_mat_pixel_IQ(acq,squeeze(rawData(:,:,an)),timeVector,focal_delay,x_grid(pixelS),z_grid(pixelS),probe_geometry,Angl,c0, w0,Fnum);
end

envApS = [];
alphaR_ApS = [];
alphaROpti = [];
tic
for i1 = 1:ItApE-31
    ApS = ItApFull(i1):ItApFull(i1+31);
    for k = 1:length(angI)
        % envApS (:,k,i1) =  abs(sum(ScatterMatrixPW(:,ApS,k),2));
        envApS_Temp (:,k) =  abs(sum(ScatterMatrixPW(:,ApS,k),2));
        
    end
    center_ApS(i1) = ceil(median(ApS));
    alphaR_ApS(i1,:) = squeeze(alphaR(:, center_ApS(i1)));
    envApS(:,:,i1)=envApS_Temp;
end

[~,maxSAIdx] = (max(envApS,[],3));

[maxEnergyVal,OptiTxIdx] = max(max(envApS,[],3),[],2);
for p = 1:length(pixelS)
    alphaROpti(p) = alphaR_ApS(maxSAIdx(p,OptiTxIdx(p)),p);
    meanP(p) = mean(ScatterMatrixPW(p,:,OptiTxIdx(p)));
    stdP(p) = std(ScatterMatrixPW(p,:,OptiTxIdx(p)));
end

mstd_R = abs(meanP./stdP);

s1(pixelS) = sind(alphaROpti(:));
c1(pixelS) = cosd(alphaROpti(:));

uPW(pixelS) = maxEnergyVal.*s1(pixelS);
vPW(pixelS) = maxEnergyVal.*c1(pixelS);

uPW(pixelS(mstd_R<threshold)) = 0;
vPW(pixelS(mstd_R<threshold)) = 0;

% uPW_interp = interp1(pixelS, uPW(pixelS), pixelFull,'linear',0);
% vPW_interp = interp1(pixelS, vPW(pixelS), pixelFull,'linear',0);
% 
% uPW(pixelFull) = uPW_interp;
% vPW(pixelFull) = vPW_interp;

uPW = reshape(uPW, size(x_grid,1),size(x_grid,2));
vPW = reshape(vPW, size(x_grid,1),size(x_grid,2));

%%
zp = 1:length(z_axis_IQ);
xp = 1:length(probe_geometry);

%figure
ax3= ax1;
colormap(ax3,'hot');
hold on,

uPWs = uPW;%smoothn(uPW, 'robust');
vPWs = vPW;% smoothn(vPW, 'robust');

%quiverC2D(x_grid(zp,xp).*100,z_grid(zp,xp).*100,uPW(zp,xp),-vPW(zp,xp),8,'LineWidth', 3,'Parent',ax2) ;
quiverC2D(x_grid(zp,xp).*1000,z_grid(zp,xp).*1000,uPWs(zp,xp),-vPWs(zp,xp),abs(beamformedDataDAS),15,'LineWidth', 3,'Parent',ax3) ;

set(gca,'TickLabelInterpreter','latex');
set(gca,'YDir','Reverse');
colormap(ax1,'gray');
xlabel(ax3,'\bf{mm}')
set(gca,'fontsize',24);
%set(gcf, 'Position',  [100, 100, 600, 600])
%pause(0.1)
hold off
toc
end




