% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: Plot1D_CurrentEstimate.m                                    %
%    Purpose: Approximates the current at the tip with the assumption of  %
%             a known, constant propagation speed defined by whether the  %
%             discharge is a streamer or leader. Outputs to the screen    %
%             relevant values and associated statistics. Creates the      %
%             following plots:                                            %
%             (1a) a highlight of the longest path and                    %
%             (1b) initiation-to-ground path for CGs,                     %
%             (2)  current estimation over the course of the discharge,   %
%             (3)  the currents w.r.t. altitude, and                      %
%             (4)  the timescale of the propagation.                      %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 3, 2024                                            %
%    Updates:     June 2025 - Introduced ability to distinguish between   %
%                             negative and positive propagation.          %
%              October 2025 - Integrated checkMagnitude.m function.       %
%             December 2025 - Integrated processing functions.            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close all
clearvars -except sims

%% Fully-automated onwards:
fprintf('\n*** Executing Plot1D_CurrentEstimate.m script. ***\n');  
if ~exist('sims','var')
    specifySimDetails;
end

% Read-in data for simulation and process the current at the tip:
[current, polarity, tau, links, rhos] = processTipCurrent(sims);

% Determine the scale:
temporalFactor    = checkMagnitude([links.timescale.min; links.timescale.max]);
currentFactor     = checkMagnitude(current.case(:));
minCurrentFactor  = checkMagnitude(current.min.val);
conversionFactor  = (10^(-9))*sims.domain.dx*sims.domain.dy*sims.domain.dz;
chargeTransported = conversionFactor.*abs(current.polarity(1:links.total));

% Output to the screen the results:
fprintf(strcat("\n\t",sims.objectType," (positive) propagation speed approximated as ",num2str(sims.vprop.pos/1000,'%.0f')," km/s."))
fprintf(strcat("\n\t",sims.objectType," (negative) propagation speed approximated as ",num2str(sims.vprop.neg/1000,'%.0f')," km/s."))
fprintf(strcat("\n\tDischarge timescale estimated to be between ",num2str(temporalFactor.Number*links.timescale.min,'%.1f')," - ",num2str(temporalFactor.Number*links.timescale.max,'%.1f')," ",temporalFactor.Unit,"s."))
fprintf(strcat("\n\tLongest branch traveled ",num2str(links.maxdistance*sims.spatialFactor.Number,'%.1f')," ",sims.spatialFactor.Unit,"m from the initiation point."))
fprintf(strcat("\n\tAverage currents (based on charge transfer):"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.avg.val,5),"\t",currentFactor.Unit,"A\t(mean, including all links,  sigma = ",num2str(currentFactor.Number*current.avg.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.partial.val,5),"\t",currentFactor.Unit,"A\t(mean, excluding final link, sigma = ",num2str(currentFactor.Number*current.partial.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.max.val,5),"\t",currentFactor.Unit,"A\t(maximum, link ",num2str(links.steps(current.max.index)),")"));
fprintf(strcat("\n\t\t",num2str(minCurrentFactor.Number*current.min.val,5),"\t",minCurrentFactor.Unit,"A\t(minimum, link ",num2str(links.steps(current.min.index)),")\n"));
if links.end(end,3)==0
    fprintf(strcat("\t\t",num2str(0.001*chargeTransported(end)./(links.pathdistance(end)/(polarity.neg(end)*sims.vprop.neg + polarity.pos(end)*sims.vprop.pos)),5),"\tkA\t(initiation-to-ground path)\n"))
end

%% Highlights longest path (or path to ground, if applicable):
[longestpath.start, longestpath.end] = string2link(links.path(links.maxindex),links.initiate);
highlightPath(longestpath,links,'r',sims)
exportgraphics(gcf,strcat(sims.pathPNGs,'/Path-Longest_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
if links.end(end,3)==0
    [groundpath.start, groundpath.end] = string2link(links.path(end),links.initiate);
    highlightPath(groundpath,links,'b',sims)
    exportgraphics(gcf,strcat(sims.pathPNGs,'/Path-Ground_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
end

%% Plot of instantaneous current values:
clf
hold on; box on;
scatter(links.steps(current.polarity>0),currentFactor.Number*current.case(current.polarity>0),50,'r.','DisplayName',strcat("($+$) Instantaneous Current: ",num2str(currentFactor.Number*mean(current.case(current.polarity>0)),3)," $\pm$ ",num2str(currentFactor.Number*std(current.case(current.polarity>0)),3)," ",currentFactor.LaTeX,"A"))
scatter(links.steps(current.polarity<0),currentFactor.Number*current.case(current.polarity<0),50,'b.','DisplayName',strcat("($-$) Instantaneous Current: ",num2str(currentFactor.Number*mean(current.case(current.polarity<0)),3)," $\pm$ ",num2str(currentFactor.Number*std(current.case(current.polarity<0)),3)," ",currentFactor.LaTeX,"A"))
yline(currentFactor.Number*current.max.val,'k:','DisplayName',strcat("Maximum Current: ",num2str(currentFactor.Number*current.max.val,3)," ",currentFactor.LaTeX,"A (link ",num2str(links.steps(current.max.index)),")"),'LineWidth',2);
if links.end(end,3)==0
    yline(currentFactor.Number*current.avg.val,'k--','DisplayName',strcat("Average Current: ",num2str(currentFactor.Number*current.partial.val,3)," $\pm$ ",num2str(currentFactor.Number*current.partial.std,3)," ",currentFactor.LaTeX,"A"),'LineWidth',2);
else
    yline(currentFactor.Number*current.avg.val,'k--','DisplayName',strcat("Average Current: ",num2str(currentFactor.Number*current.avg.val,3)," $\pm$ ",num2str(currentFactor.Number*current.avg.std,3)," ",currentFactor.LaTeX,"A"),'LineWidth',2);
end
yline(currentFactor.Number*current.min.val,'k-','DisplayName',strcat("Minimum Current: ",num2str(minCurrentFactor.Number*current.min.val,3)," ",minCurrentFactor.LaTeX,"A (link ",num2str(links.steps(current.min.index)),")"),'LineWidth',2);
xlim([0,links.steps(end)]);
legend('Interpreter','latex','FontSize',16,'location','best','box','off')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
title(strcat("Current Estimate for ",sims.objectType," Discharge on ",sims.objectName),'Interpreter','latex','FontSize',28);
xlabel('Propagation Link $l$','FontSize',20,'Interpreter','latex');
ylabel(strcat('Current $I$ (',currentFactor.LaTeX,'A)'),'FontSize',20,'Interpreter','latex');
grid off; hold off;
set(gcf,'Position',[0,0,1200,800]);
exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentEstimate_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);  

%% Plot of currents wrt end node altitudes:
clf
hold on; box on;
x.pos = currentFactor.Number*current.case(current.polarity>0);
y.pos = sims.spatialFactor.Number*links.height(current.polarity>0);
x.neg = currentFactor.Number*current.case(current.polarity<0);
y.neg = sims.spatialFactor.Number*links.height(current.polarity<0);
if y.pos(end) == 0
    scatter(x.pos(1:(end-1)),y.pos(1:(end-1)),50,'r.','DisplayName',strcat("($+$) Instantaneous Current: ",num2str(mean(x.pos(1:(end-1))),3)," $\pm$ ",num2str(std(x.pos(1:(end-1))),3)," ",currentFactor.LaTeX,"A"));
else
    scatter(x.pos,y.pos,50,'r.','DisplayName',strcat("($+$) Instantaneous Current: ",num2str(mean(x.pos),3)," $\pm$ ",num2str(std(x.pos),3)," ",currentFactor.LaTeX,"A"));
end
if y.neg(end) == 0
    scatter(x.neg(1:(end-1)),y.neg(1:(end-1)),50,'b.','DisplayName',strcat("($-$) Instantaneous Current: ",num2str(mean(x.neg(1:(end-1))),3)," $\pm$ ",num2str(std(x.neg(1:(end-1))),3)," ",currentFactor.LaTeX,"A"));
else
    scatter(x.neg,y.neg,50,'b.','DisplayName',strcat("($-$) Instantaneous Current: ",num2str(mean(x.neg),3)," $\pm$ ",num2str(std(x.neg),3)," ",currentFactor.LaTeX,"A"));
end
yline(sims.spatialFactor.Number*links.start(1,3)*sims.domain.dz,'k:','DisplayName',strcat("Initiation Height: ",num2str(sims.spatialFactor.Number*links.start(1,3)*sims.domain.dz,3)," ",sims.spatialFactor.LaTeX,"m"),'LineWidth',1);
if y.pos(end) == 0 || y.neg(end)==0
    xlim([0,currentFactor.Number*max(current.case(1:(end-1)))]);
else
    xlim([0,currentFactor.Number*max(current.case)]);
end
ylim([0 sims.spatialFactor.Number*(max(links.height)+(sims.domain.dz))]);
legend('Interpreter','latex','FontSize',16,'location','best','box','off')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
title(strcat("Current Estimate for ",sims.objectType," Discharge on ",sims.objectName),'Interpreter','latex','FontSize',28);
xlabel(strcat('Current $I$ (',currentFactor.LaTeX,'A)'),'FontSize',20,'Interpreter','latex');
ylabel(strcat('Altitude of End Node (',sims.spatialFactor.LaTeX,'m)'),'FontSize',20,'Interpreter','latex');
grid off; hold off;
set(gcf,'Position',[0,0,1200,800]);
exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentHeight_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);  

%% Plot of link-addition timescales wrt step number:
clf
hold on; box on;
scatter(links.steps,temporalFactor.Number*tau.case,25,'ko','filled','MarkerFaceAlpha',0.5,'DisplayName',strcat("Propagation Timescale: $$\Delta t_n = \frac{\Delta l_n}{v_\mathrm{prop}}$$"))% where $$v_\mathrm{prop}=$$ ",num2str(testVelocity/1000)," km/s"))
yline(mean(temporalFactor.Number*tau.case),'--k',strcat(num2str(mean(temporalFactor.Number*tau.case),'%.2f')," ",temporalFactor.LaTeX,"s"),'LineWidth',1,'DisplayName',strcat("Mean Timescale: $$\bar{\Delta t}=\frac{1}{N} \sum^N_{n=1}{\Delta t_n}$$"),'Interpreter','latex','FontSize',12)
xlim([0,links.steps(end)])
ylim([0 ceil(max(temporalFactor.Number*tau.case))])
legend('Interpreter','latex','FontSize',16,'location','best','box','off')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
title(strcat("Link Timescale (",sims.objectType,", Grid Spacings $=$ ",num2str(sims.domain.dx*sims.spatialFactor.Number)," ",sims.spatialFactor.LaTeX,"m)"),'Interpreter','latex','FontSize',24);
xlabel('Propagation Link $l$','FontSize',20,'Interpreter','latex');
ylabel(strcat('Timescale $\tau$ (',temporalFactor.LaTeX,'s)'),'FontSize',20,'Interpreter','latex');
grid off
hold off
set(gcf,'Position',[0,0,600,1000]);
exportgraphics(gcf,strcat(sims.pathPNGs,'/Timescales_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);  % Before discharge

%% Functions to isolate particular types of paths:
function highlightPath(longestpath,links,color,sims)
    dxyz = [sims.domain.dx; sims.domain.dy; sims.domain.dz]*sims.spatialFactor.Number;
    clf
    hold on
    for N = 1:1:size(links.start,1)
        plot3([(links.start(N,1)-1)*dxyz(1); (links.end(N,1)-1)*dxyz(1)],[(links.start(N,2)-1)*dxyz(2); (links.end(N,2)-1)*dxyz(2)],[(links.start(N,3)-1)*dxyz(3); (links.end(N,3)-1)*dxyz(3)],'k','LineWidth',1);
    end
    for M = 1:1:size(longestpath.end,1)
        plot3([(longestpath.start(M,1)-1)*dxyz(1); (longestpath.end(M,1)-1)*dxyz(1)],[(longestpath.start(M,2)-1)*dxyz(2); (longestpath.end(M,2)-1)*dxyz(2)],[(longestpath.start(M,3)-1)*dxyz(3); (longestpath.end(M,3)-1)*dxyz(3)],color,'LineWidth',4);
    end
    setUpAxes(sims,'xyz')
    set(gcf,'Position',[0,0,1000,1000]);
    set(gcf,'Resize','off');
end