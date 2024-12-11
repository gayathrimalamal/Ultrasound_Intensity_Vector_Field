%-- Implements the DAS RF beamforming with pixel wise addition
function beamformedData = DAS_pixel_IQ(acq,rawData,timeVector,focal_delay,x_grid,z_grid,probe_geometry,ang,c_map, w0,Fnum)
N_elements=size(rawData,2);
rows=size(x_grid,1);
columns=size(x_grid,2);
beamformedData = zeros(rows,columns);
delay_compensation=zeros(rows*columns,N_elements);

if(strcmp(acq, 'PW'))
    transmit_delay = z_grid(:)*cos(ang)+x_grid(:)*sin(ang);
else
    transmit_delay = sqrt((x_grid(:)-ang).^2+z_grid(:).^2); %ang will be the transmit center in SA transmission
end
%depth based apodization
rx_f_number = Fnum;
rx_aperture = z_grid(:)/rx_f_number;
rx_aperture_distance = abs(x_grid(:)*ones(1,N_elements)-ones(rows*columns,1)*probe_geometry.'); %abs(x_grid(:)*ones(1,N_elements)-ones(length(x_grid(:),1))*probe_geometry.');
receive_apodization = apodization(rx_aperture_distance,rx_aperture*ones(1,N_elements),'hanning');

%delay compensation
for nrx=1:N_elements
    
    receive_delay=sqrt((probe_geometry(nrx,1)-x_grid(:)).^2+z_grid(:).^2);
    delay = (transmit_delay+receive_delay)./c_map(:);
    phase_shift = exp(1i.*w0*(delay-2*z_grid(:)./c_map(:)));
    
    delay_compensation(:,nrx) = phase_shift.*interp1(timeVector+focal_delay(nrx),rawData(:,nrx),delay,'spline',0);
end
beamformedData(:)=sum(receive_apodization.*delay_compensation,2);
beamformedData(isnan(beamformedData))=0;
beamformedData=reshape(beamformedData,size(x_grid));