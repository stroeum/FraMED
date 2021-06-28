close all
clear all
clc

Q = load('CarriedCharge.dat');
step = (0:size(Q)-1)';

set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

plot(step,Q,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('Q_{cha} (C)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(Q)]) ;
box on
grid on