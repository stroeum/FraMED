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
hgexport(gcf,'ChargeTransfer.eps');
cd ../viz

%% Channel Potential
Plot1D_ChannelPotential;
cd ../Figures
hgexport(gcf,'ChannelPotential.eps');
cd ../viz

%% Dipole Moment
Plot1D_DipoleMoment;
cd ../Figures
hgexport(gcf,'DipoleMoment.eps');
cd ../viz

%% Electrostatic Energy
Plot1D_EsEnergy;
cd ../Figures
hgexport(gcf,'EsEnergy.eps');
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
hgexport(1,'InitiationRequirements.eps')
hgexport(2,'phi.eps');
hgexport(3,'E.eps');
cd ../viz

%% Field Evolution (cont.)
Plot2D_Fields3x3;
cd ../Figures
hgexport(gcf,'FieldSpaceEvolution.eps');
cd ../viz

%% Fieldlines
Plot2D_FieldLines('~');
cd ../Figures
hgexport(gcf,'FieldLines.eps');
cd ../viz

%% LMA
LMA_main;
cd ../Figures
hgexport(gcf,'LMA.eps');
cd ../viz

