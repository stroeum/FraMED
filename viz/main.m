% main file
% Plot all figures and store them as .eps
% NB: Videos should be store separately

!mkdir ../Figures
%% Clear Folder
cd ../Figures
if ~exist('PNGs', 'dir')
    mkdir('PNGs')
end
cd PNGs
!rm -rf *.png
cd ..

if ~exist('EPSs', 'dir')
    mkdir('EPSs')
end
cd EPSs
!rm -rf *.eps
cd ..
cd ../viz

%% Charge Transfer
Plot1D_ChannelCharge;
cd ../Figures/EPSs
exportgraphics(gcf,'ChargeTransfer.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'ChargeTransfer.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Channel Potential
Plot1D_ChannelPotential;
cd ../Figures/EPSs
exportgraphics(gcf,'ChannelPotential.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'ChannelPotential.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Dipole Moment
Plot1D_DipoleMoment;
cd ../Figures/EPSs
exportgraphics(gcf,'DipoleMoment.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'DipoleMoment.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Electrostatic Energy
Plot1D_EsEnergy;
cd ../Figures/EPSs
exportgraphics(gcf,'EsEnergy.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'EsEnergy.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Field Evolution
Plot1D_FieldEvolution;
cd ../Figures/EPSs
hgexport(1,'FieldTimeEvolution.eps');
hgexport(2,'phiTimeEvolution.eps');
hgexport(3,'ETimeEvolution.eps');
cd ../../viz

%% Field Evolution (cont.)
Plot1D_Fields;
cd ../Figures/EPSs
hgexport(1,'InitiationRequirements.eps');
hgexport(2,'phi.eps');
hgexport(3,'E.eps');
cd ../../viz

%% Field Evolution (cont.)
Plot2D_Fields3x3;
cd ../Figures/EPSs
exportgraphics(gcf,'FieldSpaceEvolution.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'FieldSpaceEvolution.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Fieldlines
Plot2D_FieldLines('~');
cd ../Figures/EPSs
exportgraphics(gcf,'FieldLines.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'FieldLines.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Charged cloud structure
Plot3D_CloudDistribution;
cd ../Figures/EPSs
exportgraphics(gcf,'CloudDistribution.eps','BackgroundColor','white');
cd ../PNGs
exportgraphics(gcf,'CloudDistribution.png','BackgroundColor','white','Resolution',300);
cd ../../viz

%% Lightning visualization
LightningVisual;

%% LMA
LMA_main;
cd ../Figures/EPSs
exportgraphics(gcf,'LMA.eps','BackgroundColor','white');
cd ../../viz
