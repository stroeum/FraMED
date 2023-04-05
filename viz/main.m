% main file
% Plot all figures and store them as .eps
% NB: Videos should be store separately
clc
clearvars
clf

%% Folders for Saving (User-Defined)
fprintf('*** Executing main.m script. ***\n');
prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
sims.objectName = input(prompt1,'s');
prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
sims.objectType = input(prompt2,'s');
while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
    fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
    sims.objectType = input(prompt2,'s');
end

% Ensures directories exist
sims.pathPNGs = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/PNGs');
if ~exist(sims.pathPNGs, 'dir')
    mkdir(sims.pathPNGs)
end
sims.pathEPSs = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/EPSs');
if ~exist(sims.pathEPSs, 'dir')
    mkdir(sims.pathEPSs)
end

sims.pathVideos = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/Videos');
if ~exist(sims.pathVideos, 'dir')
    mkdir(sims.pathVideos)
end

%% Charge Transfer
Plot1D_ChannelCharge;
exportgraphics(gcf,[sims.pathEPSs,'/ChargeTransfer.eps'],'BackgroundColor','white');

%% Channel Potential
Plot1D_ChannelPotential;
exportgraphics(gcf,[sims.pathEPSs,'/ChannelPotential.eps'],'BackgroundColor','white');

%% Dipole Moment
Plot1D_DipoleMoment;
exportgraphics(gcf,[sims.pathEPSs,'/DipoleMoment.eps'],'BackgroundColor','white');

%% Electrostatic Energy
Plot1D_EsEnergy;
exportgraphics(gcf,[sims.pathEPSs,'/EsEnergy.eps'],'BackgroundColor','white');

%% Field Evolution
Plot1D_FieldEvolution;
exportgraphics(figure(1),[sims.pathEPSs,'/FieldTimeEvolution.eps'],'BackgroundColor','white');
exportgraphics(figure(2),[sims.pathEPSs,'/phiTimeEvolution.eps'],'BackgroundColor','white');
exportgraphics(figure(3),[sims.pathEPSs,'/ETimeEvolution.eps'],'BackgroundColor','white');

%% Field Evolution (cont.)
Plot1D_Fields;
exportgraphics(figure(1),[sims.pathEPSs,'/InitiationRequirements.eps'],'BackgroundColor','white');
exportgraphics(figure(2),[sims.pathEPSs,'/phi.eps'],'BackgroundColor','white');
exportgraphics(figure(3),[sims.pathEPSs,'/E.eps'],'BackgroundColor','white');

%% Field Evolution (cont.)
Plot2D_Fields3x3;
exportgraphics(gcf,[sims.pathEPSs,'/FieldSpaceEvolution.eps'],'BackgroundColor','white');

%% Fieldlines
Plot2D_FieldLines('~',sims);
exportgraphics(gcf,[sims.pathEPSs,'/FieldLines.eps'],'BackgroundColor','white');

%% Charged cloud structure
Plot3D_CloudDistribution;
exportgraphics(gcf,[sims.pathEPSs,'/CloudDistribution.eps'],'BackgroundColor','white');

%% Lightning visualization
LightningVisual;

%% LMA
LMA_main;
exportgraphics(gcf,[sims.pathEPSs,'/LMA.eps'],'BackgroundColor','white');
