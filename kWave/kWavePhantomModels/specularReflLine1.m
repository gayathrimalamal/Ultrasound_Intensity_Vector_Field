%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Defines a two layer model of fat and tissue from an input imgtemplate
% Nx - grid points in x direction
% Ny - grid points in y direction
%  V1_Date 05-02-2022
%  Author: Gayathri M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phantom_model, sound_speed, density, alpha_coeff, Fnum] = specularReflLine1 (Nx, Ny, lineStartX, lineStartY, lineLen, lineThickness, tilt, p_grid)

%% SoSPhantom_TwoLayer3 template
phantom_model = zeros(Nx,Ny);
% phantom_model (imNorm>0.7) = 1;
% phantom_model (imNorm<0.5) = 2;
% phantom_model(pin==1) = 3;
Fnum = 1.5; % Receive number for beamforming
%% Define the medium properties
% Tissue 1
cp1                 = 1540;   	% Avg speed of sound in tissue [m/s]
rho1                = 1000;     % Density [kg/m^3]
alpha1              = 0.5;    % Absorption [dB/(MHz cm)]

% Tissue 2
cpRefl                = 3900;    % Speed of sound in  [m/s]
rhoRefl               = 2600;    % Density in [kg/m^3] of bone
alphaRefl             = 3.5;     % Absorption [dB/(MHz cm)]

% =========================================================================
% Creating reflector line
% =========================================================================

phantom =  makeLine(Nx, Ny, [round(lineStartX) round(lineStartY)], -tilt*pi/180, lineLen); %gamma = -tilt

for Nc=1:lineThickness
    phantom = phantom + makeLine(Nx, Ny, [round(lineStartX)+Nc  round(lineStartY)], -tilt*pi/180, lineLen);
end

%=========================================================================
%Creating background maps of random scatterers
%=========================================================================
% p_grid1 = bgTemplate;%%
% p_grid = randn([Nx, Ny]);
%  p_grid = p_grid2;%+2.*p_grid2;

background_map_mean1 = 1;
background_map_std1 = 0.01;
background_map = background_map_mean1 + background_map_std1 * p_grid;
% =========================================================================
% Defining medium speed of sound and density
% =========================================================================

sound_speed            = cp1*ones(Nx, Ny).*background_map;
sound_speed(phantom == 1) = cpRefl;
density                = rho1*ones(Nx, Ny).*background_map;
density(phantom== 1)= rhoRefl;

% if(tilt ~=0)
%     sound_speed = smoothn(sound_speed, 'robust');
%     density = smoothn(density, 'robust');
% end

alpha_coeff            = alpha1*ones(Nx, Ny);
alpha_coeff(phantom == 1) = alphaRefl;
end
