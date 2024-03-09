% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: Plot3D_CloudDistribution.m                                   %
% Purpose: Visualizes the charged cloud structure and ground with custom  %
%          colormap. Outputs a figure to the screen but does not save it  %
%          to a file automatically Updated February 20, 2024 to allow for %
%          free-space boundary condition distinction.                     %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: February 22, 2022                                           %
% Last Update: February 20, 2024                                          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Initiate
close all
clearvars -except sims

clf
cd ../results
fprintf('\n*** Executing Plot3D_CloudDistribution.m script. ***\n');

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
end 
prompt_Grounded = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
is.Grounded = input(prompt_Grounded,'s');                    
while ~strcmp(is.Grounded,'FS') && ~strcmp(is.Grounded,'G')
    fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
    is.Grounded = input(prompt_Grounded,'s');
end
promptSaveAfter = '\nWould you like to save an image of the charge layer distributions after the discharge has occurred? (Y / N)\n-->';
saveAfter = input(promptSaveAfter,'s');    
while ~strcmp(saveAfter,'Y') && ~strcmp(saveAfter,'N')
    fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
    saveAfter = input(promptSaveAfter,'s');
end

%% Load data files
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
if strcmp(saveAfter,'Y')
    prompt_fileNum = '\nWhich iteration save file would you like to read-in?\nInput example: rhoAmb3d2500.dat\n-->';
    iterationFile = input(prompt_fileNum,'s');
    while ~exist(['../results/',iterationFile],'file')
        fprintf(['\t*** File ''',iterationFile,''' does not exist in the ''results'' directory. Please try again. ***\n']);
        iterationFile = input(prompt_fileNum,'s');
    end
    rho.data = load(['../results/',iterationFile],'-ascii'); % After discharge
else
    fprintf('\tUsing charge layer data from before discharge occurs.\n');
    rho.data = load('rhoAmb.dat','-ascii'); % Before discharge
end
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

cd ../viz
rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3

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
figure(1);
set(gcf,'Position',[0,0,800,600]);
hold on;

% Calls the new function that automatically recognizes charge regions:
plottingChargeRegions('white',0.4,rho,X,Y,Z);

% Represents the neutrally charged (grounded) surface:
if strcmp(is.Grounded,'G')
    P.x = [L.x 0 0 L.x]*1e-3;
    P.y = [L.y L.y 0 0]*1e-3;
    P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
    patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
end
axis equal

% Sets bounds for the axes (comment out if clouds get cut off):
%axis([0 L.x 0 L.y 15000 55000]*1e-3)                % Cropped view
axis([0 L.x 0 L.y gnd.alt 2/2*(L.z+gnd.alt)]*1e-3)  % Full span 

%% Save figure as file:
if strcmp(saveAfter,'Y')
    exportgraphics(gcf,[sims.pathPNGs,'/CloudDistribution_After_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);  % Before discharge
else
    exportgraphics(gcf,[sims.pathPNGs,'/CloudDistribution_Before_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);   % After discharge
end