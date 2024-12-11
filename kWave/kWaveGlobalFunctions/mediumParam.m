function med = mediumParam(f)

%% DEFINE THE MEDIUM PARAMETERS FOR LUNG

%Skin layer
med.vSkin=1700;
med.rhoSkin=1150;
med.alphaSkin=0.6;

med.vSoftTissue = 1540;                      % [m/s]
med.rhoSoftTissue = 1100;                    % [kg/m^3]
med.alphaSoftTissue = 0.5;%*f;                    % [db/cm/MHz]

%Air
med.vAir = 330;
med.rhoAir = 100;
med.alphaAir= 14;%*f;

%Pleura
med.vPleura = med.vSoftTissue;
med.rhoPleura = med.rhoSoftTissue;
med.alphaPleura = med.alphaSoftTissue;

%Fluid (Water at 20 degree)
med.vFluid = 1430;
med.rhoFluid = 1000;
med.alphaFluid = 0.0022;%*f;

end
