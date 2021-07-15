close all
clearvars
clc
cd ../results/
p = load('DischargeDipoleMoment.dat')*1e-3; % convert to C/km
NbOfPoints  = size(p);
NbOfPoints  = NbOfPoints(1);
step        = (0:NbOfPoints-1)';
VectorStep  = 25;                   % plot p vector every VectorStep step.
 
figure;
set(gcf,'Units','inches','OuterPosition', [20 20 20 40]/6)
% set(gca,'Units','inches','Position', [4 3 15 11]/6)

subplot(211)
absp = (p(:,1).^2+p(:,2).^2+p(:,3).^2).^.5;
% hold on
plot(step,absp,'LineWidth',1,'LineStyle','-','color','k');
% plot(step,-p(:,3),'LineWidth',1,'LineStyle','--','color','b');
% hold off
xlabel('step','FontSize',12);
ylabel('$\|\vec{p}\|$ (C.km)','Interpreter','latex','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(absp)]) ;
box on
grid on

subplot(212)
hold on
plot(step,p(:,1),'LineWidth',1,'LineStyle','-','Color','r');
plot(step,p(:,2),'LineWidth',1,'LineStyle','-','Color','g');
plot(step,p(:,3),'LineWidth',1,'LineStyle','-','Color','b');
xlabel('step','FontSize',12);
ylabel('p_{x,y,z} (C.km)','FontSize',12);
legend('p_x','p_y','p_z','Location','East');
set(gca,'FontSize',10);
axis([0 max(step) min(min(p)) max(max(p))]);
box on
grid on
hold off

cd ../viz/