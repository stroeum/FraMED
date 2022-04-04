% main file
% Plot all figures and store them as .eps
% NB: Videos should be store separately

!mkdir ../Figures
%% Clear Folder
cd ../Figures
!rm -rf *.eps
cd ../viz

%% Charge Transfer
Plot1D_ChannelCharge;
cd ../Figures
exportgraphics(gcf,'ChargeTransfer.eps','BackgroundColor','white');
cd ../viz

%% Channel Potential
Plot1D_ChannelPotential;
cd ../Figures
exportgraphics(gcf,'ChannelPotential.eps','BackgroundColor','white');
cd ../viz

%% Dipole Moment
Plot1D_DipoleMoment;
cd ../Figures
exportgraphics(gcf,'DipoleMoment.eps','BackgroundColor','white');
cd ../viz

%% Electrostatic Energy
Plot1D_EsEnergy;
cd ../Figures
exportgraphics(gcf,'EsEnergy.eps','BackgroundColor','white');
cd ../viz

%% Field Evolution
Plot1D_FieldEvolution;
cd ../Figures
hgexport(1,'FieldTimeEvolution.eps');
hgexport(2,'phiTimeEvolution.eps');
hgexport(3,'ETimeEvolution.eps');
cd ../viz

%% Field Evolution (cont.)
Plot1D_Fields;
cd ../Figures
hgexport(1,'InitiationRequirements.eps');
hgexport(2,'phi.eps');
hgexport(3,'E.eps');
cd ../viz

%% Field Evolution (cont.)
Plot2D_Fields3x3;
cd ../Figures
exportgraphics(gcf,'FieldSpaceEvolution.eps','BackgroundColor','white');
cd ../viz

%% Fieldlines
Plot2D_FieldLines('~');
cd ../Figures
exportgraphics(gcf,'FieldLines.eps','BackgroundColor','white');
cd ../viz

%% LMA
LMA_main;
cd ../Figures
exportgraphics(gcf,'LMA.eps','BackgroundColor','white');
cd ../viz

%% Charged cloud structure
Plot3D_CloudDistribution;
cd ../Figures
exportgraphics(gcf,'CloudDistribution.eps','BackgroundColor','white');
cd ../viz

