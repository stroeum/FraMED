% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: Plot3D_ChargeLayerDefinitions.m                              %
% Purpose: Visualizes the charged cloud structure and ground with custom  %
%          colormap. Labels the charge layers with their FraMED input.    %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: April 29, 2022                                              %
% Last Update: N/A                                                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Initiate
close all
clearvars -except sims

clf
cd ../results

if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
    prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
    sims.objectType = input(prompt2,'s');
    while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
        fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
        sims.objectType = input(prompt2,'s');
    end

    % Settings to ensure proper directory referencing:
    sims.pathPNGs = ['../Figures/',sims.objectName,'/',sims.objectType,'/PNGs'];
    if ~exist(sims.pathPNGs,'dir')
        mkdir(sims.pathPNGs);
    end
    sims.pathVideos = ['../Figures/',sims.objectName,'/',sims.objectType,'/Videos'];
    if ~exist(sims.pathVideos,'dir')
        mkdir(sims.pathVideos);
    end

    % Specifies the boundary conditions for the simulation:
    prompt_BCtype = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
    sims.BCtype = input(prompt_BCtype,'s');                    
    while ~strcmp(sims.BCtype,'FS') && ~strcmp(sims.BCtype,'G')
        fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
        sims.BCtype = input(prompt_BCtype,'s');
    end
end 

% Quantities of charge layer 1:
prompt_Q1 = '\nWhat is the charge (in C) within the lower charge layer?\n-->';
Q1 = input(prompt_Q1);                
prompt_R1 = '\nWhat is the disk radius (in km) of the lower charge layer?\n-->';
R1 = input(prompt_R1);  
prompt_h1 = '\nWhat is the disk height (in km) of the lower charge layer?\n-->';
h1 = input(prompt_h1);  
prompt_center1 = '\nWhere is the center of the lower disk height (in km) of the lower charge layer?\nInput format: [ (x-coordinate); (y-coordinate); (z-coordinate)]\n-->';
center1 = input(prompt_center1);

% Quantities of charge layer 2:
prompt_Q2 = '\nWhat is the charge (in C) within the lower charge layer?\n-->';
Q2 = input(prompt_Q2);               
prompt_R2 = '\nWhat is the disk radius (in km) of the lower charge layer?\n-->';
R2 = input(prompt_R2);  
prompt_h2 = '\nWhat is the disk height (in km) of the lower charge layer?\n-->';
h2 = input(prompt_h2);  
prompt_center2 = '\nWhere is the center of the lower disk height (in km) of the lower charge layer?\nInput format: [ (x-coordinate); (y-coordinate); (z-coordinate)]\n-->';
center2 = input(prompt_center2);

%% Load data files
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
rho.data = load('rhoAmb.dat','-ascii');
gnd.alt  = load('z_gnd.dat', '-ascii');

%% Derive main parameters
N.x = Nxyz(1);
N.y = Nxyz(2);
N.z = Nxyz(3);

d.x = dxyz(1);           % _m
d.y = dxyz(2);           % _m
d.z = dxyz(3);           % _m

L.x = (N.x-1)*d.x;         % _m
L.y = (N.y-1)*d.y;         % _m
L.z = (N.z-1)*d.z;         % _m

rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3

clear Nxyz
clear dxyz
clear InitPoint
cd ../viz

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

% Calls the new function that automatically recognizes charge regions:
plottingLayerDefs('white',0.4,rho,X,Y,Z,R1,h1,Q1,center1,R2,h2,Q2,center2);

% Represents the neutrally charged (grounded) surface:
P.x = [L.x 0 0 L.x]*1e-3;
P.y = [L.y L.y 0 0]*1e-3;
P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'DisplayName','Ground');
title(['Charge Layer Distribution (',sims.objectName,')'],'Interpreter','latex','FontSize',28);
    
exportgraphics(gcf,[sims.pathPNGs,'/ChargeLayerDefs_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);
% view([90,0]) 
% camlight; lighting gouraud

