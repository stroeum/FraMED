close all
clear all
clc
cd ../results/
% 1 pica = 12 points = 1/6 inch
%set(gca,'Units','inches','Position',[2 2 20 15]/6,'OuterPosition', [0 0 30 25]/6)
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

phio = load('ChannelPotentials.dat')*1e-6;
step = (0:size(phio)-1)';

plot(step,phio,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('\phi_0 (MV)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(phio) max(phio)]) ;
box on
grid on
cd ../viz/