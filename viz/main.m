% main file
% Plot all figures and store them as .eps
% NB: Videos should be store separately
clc
close all
clearvars

%% Folders for Saving (User-Defined)
sims = specifySimDetails();

%% Runtime Results
Plot1D_RuntimeResults;
exportgraphics(gcf,strcat(sims.pathEPSs,'/RuntimeResults.eps'),'BackgroundColor','white');

%% Charge Transfer
Plot1D_ChannelCharge;
exportgraphics(gcf,strcat(sims.pathEPSs,'/ChargeTransfer.eps'),'BackgroundColor','white');

%% Channel Potential
Plot1D_ChannelPotential;
exportgraphics(gcf,strcat(sims.pathEPSs,'/ChannelPotential.eps'),'BackgroundColor','white');

%% Current Estimate
Plot1D_CurrentEstimate;
exportgraphics(gcf,strcat(sims.pathEPSs,'/CurrentEstimate.eps'),'BackgroundColor','white');

%% Dipole Moment
Plot1D_DipoleMoment;
exportgraphics(gcf,strcat(sims.pathEPSs,'/DipoleMoment.eps'),'BackgroundColor','white');

%% Electrostatic Energy
Plot1D_EsEnergy;
exportgraphics(gcf,strcat(sims.pathEPSs,'/EsEnergy.eps'),'BackgroundColor','white');

%% Field Evolution
Plot1D_FieldEvolution;
exportgraphics(figure(1),strcat(sims.pathEPSs,'/FieldTimeEvolution.eps'),'BackgroundColor','white');
exportgraphics(figure(2),strcat(sims.pathEPSs,'/phiTimeEvolution.eps'),'BackgroundColor','white');
exportgraphics(figure(3),strcat(sims.pathEPSs,'/ETimeEvolution.eps'),'BackgroundColor','white');

%% Field Evolution (cont.)
Plot1D_Fields;
exportgraphics(figure(1),strcat(sims.pathEPSs,'/InitiationRequirements.eps'),'BackgroundColor','white');
exportgraphics(figure(2),strcat(sims.pathEPSs,'/phi.eps'),'BackgroundColor','white');
exportgraphics(figure(3),strcat(sims.pathEPSs,'/E.eps'),'BackgroundColor','white');

%% Field Evolution (cont.)
Plot2D_Fields3x3;
exportgraphics(gcf,strcat(sims.pathEPSs,'/FieldSpaceEvolution.eps'),'BackgroundColor','white');

%% Fieldlines
Plot2D_FieldLines('~',sims);
exportgraphics(gcf,strcat(sims.pathEPSs,'/FieldLines.eps'),'BackgroundColor','white');

%% Charged cloud structure
Plot3D_CloudDistribution;
exportgraphics(gcf,strcat(sims.pathEPSs,'/CloudDistribution.eps'),'BackgroundColor','white');

%% Lightning visualization
LightningVisual;

%% LMA
LMA_main;
exportgraphics(gcf,strcat(sims.pathEPSs,'/LMA.eps'),'BackgroundColor','white');