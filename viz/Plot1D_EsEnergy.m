close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs')
    specifySimDetails;
end 

cd ../results/
Es = load('EsEnergy.dat');
if isempty(Es) || size(Es,1)-1 == 0
    fprintf('\n*** Plot1D_EsEnergy.m cannot be executed with current EsEnergy.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_EsEnergy.m script. ***\n');
end
cd ../viz
%DischargeDipoleMoment
step = (0:size(Es,1)-1)';

set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
energyFactor = checkMagnitude(Es);
plot(step,energyFactor.Number*Es,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12,'interpreter','latex');
ylabel(strcat('$\epsilon_{es}$ (',energyFactor.LaTeX,'J)'),'FontSize',12,'interpreter','latex');
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis([0 max(step) min(0) energyFactor.Number*max(Es)]) ;
box on
grid on
exportgraphics(gcf,strcat(sims.pathPNGs,'/EsEnergy_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
