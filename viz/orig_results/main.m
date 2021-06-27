% main file
% Plot all figures and store them as .eps
% NB: Videos should be store separately

%% Clear Folder
cd Figures
!rm -rf *.eps
cd ..

%% Charge Transfer
Plot1D_ChannelCharge;
cd Figures
hgexport(gcf,'ChargeTransfer.eps');
cd ..

%% Channel Potential
Plot1D_ChannelPotential;
cd Figures
hgexport(gcf,'ChannelPotential.eps');
cd ..

%% Dipole Moment
Plot1D_DipoleMoment;
cd Figures
hgexport(gcf,'DipoleMoment.eps');
cd ..

%% Electrostatic Energy
Plot1D_EsEnergy;
cd Figures
hgexport(gcf,'EsEnergy.eps');
cd ..

%% Field Evolution
Plot1D_FieldEvolution;
cd Figures
hgexport(1,'FieldTimeEvolution.eps');
hgexport(2,'phiTimeEvolution.eps');
hgexport(3,'ETimeEvolution.eps');
cd ..

%% Field Evolution (cont.)
Plot1D_Fields;
cd Figures
hgexport(2,'phi.eps');
hgexport(3,'E.eps');
cd ..

%% Field Evolution (cont.)
Plot2D_Fields3x3;
cd Figures
hgexport(gcf,'FieldSpaceEvolution.eps');
cd ..

%% Fieldlines
Plot2D_FieldLines;
cd Figures
hgexport(gcf,'FieldLines.eps');
cd ..

%% LMA
Plot2D_LMA_BETA;
cd Figures
hgexport(gcf,'LMA.eps');
cd ..

