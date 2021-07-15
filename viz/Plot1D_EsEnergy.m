close all
clearvars
clc

cd ../results/
Es = load('EsEnergy.dat')*1e-9;
step = (0:size(Es)-1)';

set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

plot(step,Es,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('\epsilon_{es} (GJ)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(Es)]) ;
box on
grid on
cd ../viz/