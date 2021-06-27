close all
clear all
clc

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

% figure;
% color = colormap(jet(NbOfPoints));
% hold on
% for ii=1:NbOfPoints-1
%     % plots ending point of each link %
%     plot3(...
%         [step(ii),step(ii+1)],...
%         [0 0],...
%         [absp(ii),absp(ii+1)],...
%         'Color',color(ii,:), 'LineWidth', 2 ...
%         );
%     if(rem(ii,VectorStep) == 0)
%         quiver3(step(ii),0,absp(ii),p(ii,1),p(ii,2),p(ii,3),0,'Color',color(ii,:), 'LineWidth', 2);
%     end
% % axis([0 max(step) min(p(:,2)) max(p(:,2)) min(p(:,3)) max(absp+p(:,3))]);
% axis([0 max(step)...
%     -max(abs(p(:,2))) max(abs(p(:,2)))...
%     min(absp+p(:,3)) max(max(absp),max(absp+p(:,3)))]);
% end
% hold off;
% xlabel('step','FontSize',16);
% zlabel('$\|\vec{p}\|$ (C.km)','Interpreter','latex','FontSize',16);
% view([.5 -1 .5])
% box off
% grid on
% quiver3([0 1 2],[0 0 0],[0 1 2], [p(1,1) p(2,1) p(3,1)],[p(1,2) p(2,2) p(3,2)],[p(1,3) p(2,3) p(3,3)], 1)
% quiver3(step,0*step,absp,  p(:,1),p(:,2),p(:,3),2)
% quiver3(N(:),0*N(:),0*N(:), p(:,1),p(:,2),p(:,3),.5)
