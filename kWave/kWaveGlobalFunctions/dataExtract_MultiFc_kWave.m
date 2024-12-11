function  dataset = dataExtract_MultiFc_kWave(filepath, experiment, dataset_name, dataset_owner, acq, Trans, tx_steering_delay, probe_geometry, RcvData, kgrid, medium, txFreqN)

if ~exist([filepath,'\kWaveDatasets\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name], 'dir')
    mkdir([filepath,'\kWaveDatasets\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name]);
end

fsavepath = [filepath,'\kWaveDatasets\', dataset_owner,'\',acq,'\',experiment,'\',dataset_name];
dataset.probe_geometry = probe_geometry;
dataset.Fs = 1/kgrid.dt;

dataset.Trans = Trans;
dataset.probe_name = Trans.name; %Probe name
dataset.channels = Trans.tx_N_elements;
dataset.c0 = Trans.tx_c0; %Speed of sound
dataset.Fc = Trans.tx_Fc(txFreqN)/1e6;
dataset.medium = medium;

%% Added for downsampling the k-wave data
freqDown = round(dataset.Fs/(5*dataset.Fc*1e6)); %The data in k_wave is heavily oversampled. To reduce the sampling rate to 5 time the transmit frequency.
dataset.Fs = (dataset.Fs*1e-6/freqDown);

%%
if(strcmp(acq, 'SA'))
    dataset.angles = 0;
    dataset.rawData = RcvData;
    save([fsavepath,'\','dataset','_',acq,num2str(length(Trans.tx_active_SA)),'_',dataset_name,'_',num2str(Trans.tx_Fc(txFreqN)/10^6),'MHz','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.mat'], 'dataset','-v7.3');
end

if(strcmp(acq, 'PW'))
    angles = Trans.tx_pw_steering_angles;
    dataset.angles = angles.*pi/180;
    
    for txAng = 1:length(angles)
        init_delay = (tx_steering_delay(:, txAng).*kgrid.dt);
        maxDelay = (max(init_delay)/2)*(1/kgrid.dt);
        
        if(maxDelay==0)
            maxDelay=1;
        end
        
        rawData1 = squeeze(RcvData(round(maxDelay):end,:,txAng));
        rawData(1:size(rawData1,1),:,txAng) = rawData1;
        %         dataset.rawData(1:size(rawData1,1),:,txAng) = resample(rawData1,1,freqDown);
    end
    for txAng = 1:length(angles)
        rawData_res(:,:,txAng) = resample(rawData(:,:,txAng),1,freqDown);
    end
    dataset.rawData = single(rawData_res);
    save([fsavepath,'\','dataset','_',acq, num2str(Trans.na),'_',dataset_name,'_',num2str(Trans.tx_Fc(txFreqN)/10^6),'MHz','_',datestr(now, 'dd_mmm_yyyy_HH_MM_SS_AM'),'.mat'], 'dataset','-v7.3');
end

