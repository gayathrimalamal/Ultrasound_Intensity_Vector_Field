%%
%-- Implements the DAS beamforming with columnwise addition
function beamformedData=DAS_column_PW(rawData,timeVector,x_grid,z_grid,probe_geometry,ang,c_map,Fnum)

sample_factor=1;

rawData = resample(double(rawData),sample_factor,1);
N_elements = size(rawData,2);
c = 1540;
rows=size(x_grid,1);
columns=size(x_grid,2);
beamformedData = zeros(rows,columns);

for data_ang=1:length(ang)
    
    beamformedDataTemp = zeros(rows,columns);
    delay_compensation=zeros(rows,N_elements);
   
    for Nc=1:columns
           %  Nc
        transmit_delay = z_grid(:,Nc)*cos(ang(data_ang))+x_grid(:,Nc)*sin(ang(data_ang));
        
%         transmit_delay = sqrt((z_grid(:,Nc)*cos(ang(data_ang))).^2+(x_grid(:,Nc)*sin(ang(data_ang))).^2);
        %depth based apodization
        rx_f_number = Fnum;
        rx_aperture = z_grid(:,Nc)/rx_f_number;
        rx_aperture_distance = abs(x_grid(:,Nc)*ones(1,N_elements)-ones(length(x_grid(:,Nc)),1)*probe_geometry.');
        receive_apodization = apodization(rx_aperture_distance,rx_aperture*ones(1,N_elements),'tukey25');
               
        %delay compnesation
        for nrx=1:N_elements
            
            receive_delay=sqrt((probe_geometry(nrx,1)-x_grid(:,Nc)).^2+z_grid(:,Nc).^2);
            delay = (transmit_delay+receive_delay)./c;
            %             delayIdx(:,nrx,Nc)=delay/median(gradient(timeVector));
            
            %delay_compensation(:,nrx) = interp1(timeVector+focal_delay(nrx),rawData(:,nrx,data_ang),delay,'spline',0);
                         delay_compensation(:,nrx) = interp1(timeVector,rawData(:,nrx,data_ang),delay,'spline',0);
        end
        beamformedDataTemp(:,Nc)=sum(receive_apodization.*delay_compensation,2);
    end
    beamformedData=beamformedData+(1/length(ang))*beamformedDataTemp;
    
end
beamformedData(isnan(beamformedData))=0;

beamformedData=reshape(beamformedData,size(x_grid));


 