% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   File Name: Plot3D_ChargeLayerDefinitions.m                            %
%     Purpose: Visualizes the charged cloud structure and ground with     %
%              custom colormap. Labels the charge layers based on the     %
%              inputs for FraMED.                                         %
%      Author: Annelisa Esparza                                           %
%     Contact: annelisa.esparza@my.erau.edu                               %
%  Added Date: April 29, 2022                                             %
% Last Update: February 2025 - Updated to match the options available for %
%                              the createRedBlueColorMap function.        %                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Initiate
close all
clearvars -except sims

clf

if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    sims = specifySimDetails();
end

%% Load data files
cd ../results
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
rho.data = load('rhoAmb.dat','-ascii');
gnd.alt  = load('z_gnd.dat', '-ascii');
load('ChargeLayers.dat',     '-ascii');
cd ../viz

Q1      = ChargeLayers(1,1);
R1      = ChargeLayers(1,5)/1000;
h1      = ChargeLayers(1,7)/1000;
center1 = [ChargeLayers(1,2); ChargeLayers(1,3); ChargeLayers(1,4)]/1000;

Q2      = ChargeLayers(2,1);
R2      = ChargeLayers(2,5)/1000;
h2      = ChargeLayers(2,7)/1000;
center2 = [ChargeLayers(2,2); ChargeLayers(2,3); ChargeLayers(2,4)]/1000;

%% Derive main parameters
N.x = Nxyz(1);        N.y = Nxyz(2);        N.z = Nxyz(3);
d.x = dxyz(1);        d.y = dxyz(2);        d.z = dxyz(3);      % in meters
L.x = (N.x-1)*d.x;    L.y = (N.y-1)*d.y;    L.z = (N.z-1)*d.z;  % in meters

rho.data = ConvertTo3d(rho.data,Nxyz); % _nC/_m^3

clear Nxyz
clear dxyz
clear InitPoint

x = ((0:(N.x-1))*d.x)*1e-3;
y = ((0:(N.y-1))*d.y)*1e-3;
z = ((0:(N.z-1))*d.z + gnd.alt)*1e-3;

[X,Y,Z] = meshgrid(x,y,z);
rho.max = max(max(max(rho.data)));
rho.min = min(min(min(rho.data)));

%% Plotting figure:
% Map ColorScale
gnd.color = [.75 .75 .75]; % light gray for neutral charges
% Draw the tree
figure(1);
set(gcf,'Position',[0,0,1000,1000]);
hold on;

% Sets bounds for the axes (comment out if clouds get cut off):
%axis([L.x*1/5 L.x*4/5 L.y*1/5 L.y*4/5 gnd.alt 2/2*(L.z+gnd.alt)]*1e-3) % Slight crop
axis([0 L.x 0 L.y gnd.alt 2/2*(L.z+gnd.alt)]*1e-3)                     % Full span 

% Represents the neutrally charged (grounded) surface:
if strcmp(sims.BCtype,'G')
    P.x = [L.x 0 0 L.x]*1e-3;
    P.y = [L.y L.y 0 0]*1e-3;
    P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
    patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'DisplayName',strcat("Ground: $z$ = ",num2str(gnd.alt*1e-3)," km"));
end
title(['Charge Layer Distribution (',sims.objectType,' on ',sims.objectName,')'],'Interpreter','latex','FontSize',28);
   
% Calls the new function that automatically recognizes charge regions:
plottingLayerDefs('scalar',0.4,rho,X,Y,Z,R1,h1,Q1,center1,R2,h2,Q2,center2);
exportgraphics(gcf,[sims.pathPNGs,'/ChargeLayerDefs_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);

