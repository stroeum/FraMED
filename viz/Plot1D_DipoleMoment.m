close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    specifySimDetails;
end

cd ../results/
p = load('DischargeDipoleMoment.dat');
if isempty(p) || (size(p,1)-1) == 0
    fprintf('\n*** Plot1D_DipoleMoment.m cannot be executed with current DischargeDipoleMoment.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_DipoleMoment.m script. ***\n');
end
cd ../viz/

momentFactor = checkMagnitude(p);
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
plot(step,momentFactor.Number*absp,'LineWidth',1,'LineStyle','-','color','k');
% plot(step,-momentFactor.Number*p(:,3),'LineWidth',1,'LineStyle','--','color','b');
% hold off
xlabel('step','Interpreter','latex','FontSize',12);
ylabel(strcat('$\|\vec{p}\|$ (',momentFactor.LaTeX,'C$\cdot$m)'),'Interpreter','latex','FontSize',12);
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis([0 max(step) min(0) momentFactor.Number*max(absp)]) ;
box on
grid on

subplot(212)
hold on
plot(step,momentFactor.Number*p(:,1),'LineWidth',1,'LineStyle','-','Color','r');
plot(step,momentFactor.Number*p(:,2),'LineWidth',1,'LineStyle','-','Color','g');
plot(step,momentFactor.Number*p(:,3),'LineWidth',1,'LineStyle','-','Color','b');
xlabel('step','Interpreter','latex','FontSize',12);
ylabel(strcat('$p_{x,y,z}$ (',momentFactor.LaTeX,'C$\cdot$m)'),'Interpreter','latex','FontSize',12);
legend('$p_x$','$p_y$','$p_z$','Location','East','interpreter','latex');
set(gca,'FontSize',10,'TickLabelInterpreter','latex')
axis([0 max(step) momentFactor.Number*min(min(p)) momentFactor.Number*max(max(p))]);
box on
grid on
hold off
exportgraphics(gcf,strcat(sims.pathPNGs,'/DipoleMoment_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
