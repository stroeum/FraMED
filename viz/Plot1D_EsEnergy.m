close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs')
    sims = specifySimDetails();
end 

cd ../results/
Es = load('EsEnergy.dat')*1e-9;
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

plot(step,Es,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('\epsilon_{es} (GJ)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(Es)]) ;
box on
grid on
exportgraphics(gcf,strcat(sims.pathPNGs,'/EsEnergy_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);


cd ../viz/