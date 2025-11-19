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
%    Updates:    June 2025 - Introduced ability to distinguish between    %
%                            negative and positive propagation.           %
%             October 2025 - Integrated checkMagnitude.m function.        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close all
clearvars -except sims

fprintf('\n*** Executing Plot1D_CurrentEstimate.m script. ***\n');  
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    specifySimDetails;
end

%% Load data from results folder:
cd ../results/
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
load('EstablishedLinks.dat', '-ascii');
TransportedRho = load('TransportedRhoPos.dat','-ascii');
TransportedRhoCheck = load('TransportedRhoNeg.dat','-ascii');
current.polarity = load('TransportedRhoEnd.dat','-ascii');

%% Considers the case for a CG strike depositing charge to the ground:
if current.polarity(end)==0 && EstablishedLinks(end,6)==0
    current.polarity(end)=-(TransportedRho(end)+TransportedRhoCheck(end));
% Considers a now-resolved bug that recorded an additional link for ICs:
elseif abs(current.polarity(end)/current.polarity(end-1))<=(1e-14) && size(current.polarity,1)==(size(EstablishedLinks,1)+1)
    TransportedRho(end)=[];
    TransportedRhoCheck(end)=[];
    current.polarity(end)=[];
end
cd ../viz/

% Assign important values regarding the links:
links.start    = EstablishedLinks(:,1:3);
links.end      = EstablishedLinks(:,4:6);
links.distance = EstablishedLinks(:,7);
links.initiate = links.start(1,:);
links.steps    = 1:1:size(EstablishedLinks,1);
links.height   = EstablishedLinks(:,6)*dxyz(3);
spatialFactor  = checkMagnitude(links.height);

%% Calculate current:
% Convert charge density (nC/m3) into charge (C):
conversionFactor  = (10^(-9))*dxyz(1)*dxyz(2)*dxyz(3);
chargeTransported = conversionFactor.*abs(current.polarity(1:links.steps(end)));

% Isolate polarities:
polarity.pos = current.polarity;
polarity.pos(current.polarity>0) = 1;
polarity.pos(current.polarity<0) = 0;

polarity.neg = current.polarity;
polarity.neg(current.polarity>0) = 0;
polarity.neg(current.polarity<0) = 1;

% Declared constants:
vprop.leader.pos   = 4.4*(10^5); % propagation speed for positive leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
vprop.leader.neg   = 4.4*(10^5); % propagation speed for negative leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
vprop.streamer.pos = 3.0*(10^7); % propagation speed for positive streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9
vprop.streamer.neg = 3.0*(10^7); % propagation speed for negative streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9

% Calculates instantaneous timescale on a per-link basis:
tau.pos  = polarity.pos.*EstablishedLinks(:,7)./sims.vprop.pos;
tau.neg  = polarity.neg.*EstablishedLinks(:,7)./sims.vprop.neg;
tau.case = tau.pos + tau.neg;

% Calculates the current (units of Amps):
current.pos  = polarity.pos.*chargeTransported./(EstablishedLinks(:,7)./sims.vprop.pos);
current.neg  = polarity.neg.*chargeTransported./(EstablishedLinks(:,7)./sims.vprop.neg);
current.case = current.pos + current.neg;

% Reduces sets of polarity-isolated current values to nonzero values:
current.pos = nonzeros(current.pos);
current.neg = nonzeros(current.neg);

%% Identifies total branch lengths to predict lower limit of timescale:
% Determines various path branches and distances along each branch:
statusBar = waitbar(0,strcat("(",num2str(0,'%.0f'),"%) Processing ",num2str(0)," out of ",num2str(size(EstablishedLinks,1))," links."));
links.total = size(EstablishedLinks,1);
for N = 1:1:links.total
    waitbar((N/links.total),statusBar,strcat("(",num2str(100*(N-1)/links.total,'%.0f'),"%) Processing ",num2str(N)," out of ",num2str(size(EstablishedLinks,1))," links."));
    if N == 1
        links.startnode = 1;
        links.endnode = 2;
    else
        links.endnode = N+1;
    end
    
    % If the starting node is the initiation point:
    if links.start(N,:) == links.initiate
        % Classify path for this particular link:
        links.path(N) = string(characterizePropagation(links,N));
        % Total charge along this path/branch:
        links.pathcharge(N) = chargeTransported(N);
        % Total distance along this path/branch:
        links.pathdistance(N) = links.distance(N);
        % Track total distance of all links:
        if N == 1
            links.totaldistance(N) = links.distance(N);
        else
            links.totaldistance(N) = links.totaldistance(N-1)+links.distance(N);
        end
    % Otherwise, if the start node is not at the initiation point...
    else
        % Determine which path/branch this link is an extension of:
        match = find(ismember(links.end(1:(N-1),:),links.start(N,:),'row'));
        % Classify path for this particular link:
        links.path(N) = string(strcat(string(links.path(match)),string(characterizePropagation(links,N))));
        reshape(links.path,[N,1]);
        % Total distance along this path/branch:
        links.pathdistance(N) = links.pathdistance(match) + links.distance(N);
        % Total charge carried along this path/branch:
        links.pathcharge(N) = links.pathcharge(match) + chargeTransported(N);
        % Track total distance of all links:
        links.totaldistance(N) = links.totaldistance(N-1)+links.distance(N);
    end
end
close(statusBar);
[links.maxdistance, maxindex] = max(links.pathdistance);
links.timescale.max = sum(tau.case);
links.timescale.min = links.maxdistance/(polarity.neg(maxindex)*sims.vprop.neg + polarity.pos(maxindex)*sims.vprop.pos);

% Highlights longest path (or path to ground, if applicable):
[longestpath.start, longestpath.end] = string2link(links.path(maxindex),links.initiate);
highlightPath(longestpath,links,Nxyz,dxyz,'r',spatialFactor)
exportgraphics(gcf,strcat(sims.pathPNGs,'/Path-Longest_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
if links.end(end,3)==0
    [groundpath.start, groundpath.end] = string2link(links.path(end),links.initiate);
    highlightPath(groundpath,links,Nxyz,dxyz,'b',spatialFactor)
    exportgraphics(gcf,strcat(sims.pathPNGs,'/Path-Ground_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
end

%% Screen print-out of findings:
[current.min.val, current.min.index] = min(current.case);
[current.max.val, current.max.index] = max(current.case);
current.avg.val = mean(current.case); 
current.avg.std = std(current.case);
current.avg.pos = mean(current.pos);
current.avg.neg = mean(current.neg);
current.partial.val = mean(current.case(1:(end-1)));
current.partial.pos = mean(nonzeros(current.case(1:(end-1)).*polarity.pos(1:(end-1))));
current.partial.neg = mean(nonzeros(current.case(1:(end-1)).*polarity.neg(1:(end-1))));
current.partial.std = std(current.case(1:(end-1)));

% Determine the scale:
temporalFactor   = checkMagnitude([links.timescale.min; links.timescale.max]);
currentFactor    = checkMagnitude(current.case(:));
minCurrentFactor = checkMagnitude(current.min.val);

% Output to the screen the results:
fprintf(strcat("\n\t",sims.objectType," (positive) propagation speed approximated as ",num2str(sims.vprop.pos/1000,'%.0f')," km/s."))
fprintf(strcat("\n\t",sims.objectType," (negative) propagation speed approximated as ",num2str(sims.vprop.neg/1000,'%.0f')," km/s."))
fprintf(strcat("\n\tDischarge timescale estimated to be between ",num2str(temporalFactor.Number*links.timescale.min,'%.1f')," - ",num2str(temporalFactor.Number*links.timescale.max,'%.1f')," ",temporalFactor.Unit,"s."))
fprintf(strcat("\n\tLongest branch traveled ",num2str(links.maxdistance*spatialFactor.Number,'%.1f')," ",spatialFactor.Unit,"m from the initiation point."))
fprintf(strcat("\n\tAverage currents (based on charge transfer):"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.avg.val,5),"\t",currentFactor.Unit,"A\t(mean, including all links,  sigma = ",num2str(currentFactor.Number*current.avg.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.partial.val,5),"\t",currentFactor.Unit,"A\t(mean, excluding final link, sigma = ",num2str(currentFactor.Number*current.partial.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.max.val,5),"\t",currentFactor.Unit,"A\t(maximum, link ",num2str(links.steps(current.max.index)),")"));
fprintf(strcat("\n\t\t",num2str(minCurrentFactor.Number*current.min.val,5),"\t",minCurrentFactor.Unit,"A\t(minimum, link ",num2str(links.steps(current.min.index)),")\n"));
if links.end(end,3)==0
    fprintf(strcat("\t\t",num2str(0.001*chargeTransported(end)./(links.pathdistance(end)/(polarity.neg(end)*sims.vprop.neg + polarity.pos(end)*sims.vprop.pos)),5),"\tkA\t(initiation-to-ground path)\n"))
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
y.pos = spatialFactor.Number*links.height(current.polarity>0);
x.neg = currentFactor.Number*current.case(current.polarity<0);
y.neg = spatialFactor.Number*links.height(current.polarity<0);
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
yline(spatialFactor.Number*links.start(1,3)*dxyz(3),'k:','DisplayName',strcat("Initiation Height: ",num2str(spatialFactor.Number*links.start(1,3)*dxyz(3),3)," ",spatialFactor.LaTeX,"m"),'LineWidth',1);
if y.pos(end) == 0 || y.neg(end)==0
    xlim([0,currentFactor.Number*max(current.case(1:(end-1)))]);
else
    xlim([0,currentFactor.Number*max(current.case)]);
end
ylim([0 spatialFactor.Number*(max(links.height)+(dxyz(3)))]);
legend('Interpreter','latex','FontSize',16,'location','best','box','off')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
title(strcat("Current Estimate for ",sims.objectType," Discharge on ",sims.objectName),'Interpreter','latex','FontSize',28);
xlabel(strcat('Current $I$ (',currentFactor.LaTeX,'A)'),'FontSize',20,'Interpreter','latex');
ylabel(strcat('Altitude of End Node (',spatialFactor.LaTeX,'m)'),'FontSize',20,'Interpreter','latex');
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
title(strcat("Link Timescale (",sims.objectType,", Grid Spacings $=$ ",num2str(dxyz(1)*spatialFactor.Number)," ",spatialFactor.LaTeX,"m)"),'Interpreter','latex','FontSize',24);
xlabel('Propagation Link $l$','FontSize',20,'Interpreter','latex');
ylabel(strcat('Timescale $\tau$ (',temporalFactor.LaTeX,'s)'),'FontSize',20,'Interpreter','latex');
grid off
hold off
set(gcf,'Position',[0,0,600,1000]);
exportgraphics(gcf,strcat(sims.pathPNGs,'/Timescales_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);  % Before discharge


%% Functions to isolate particular types of paths:
function letter = characterizePropagation(links, N)
    % Letters A - I:
    if links.end(N,1) == (links.start(N,1)+1)
        % Letters A - C:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'A';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'B';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'C';
            end
        % Letters D - F:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'D';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'E';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'F';
            end
        % Letters G - I:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'G';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'H';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'I';
            end
        end
    % Letters J - Q:
    elseif links.end(N,1) == links.start(N,1)
        % Letters J - L:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'J';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'K';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'L';
            end
        % Letters M - N:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'M';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'N';
            end
        % Letters O - Q:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'O';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'P';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'Q';
            end
        end
    % Letters R - Z:
    elseif links.end(N,1) == (links.start(N,1)-1)
        % Letters R - T:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'R';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'S';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'T';
            end
        % Letters U - W:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'U';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'V';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'W';
            end
        % Letters X - Z:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'X';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'Y';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'Z';
            end
        end
    end
end

function [starts, ends] = string2link(path, initiation)
    currentpath = convertStringsToChars(path);
    starts = zeros([length(currentpath),3]);
    ends = starts;
    starts(1,:) = initiation;
    for N = 1:1:length(currentpath)
        switch currentpath(N)
            case 'A'
                direction = [1 -1 1];
            case 'B'
                direction = [1 0 1];
            case 'C'
                direction = [1 1 1];
            case 'D'
                direction = [1 -1 0];
            case 'E'
                direction = [1 0 0];
            case 'F'
                direction = [1 1 0];
            case 'G'
                direction = [1 -1 -1];
            case 'H'
                direction = [1 0 -1];
            case 'I'
                direction = [1 1 -1];
            case 'J'
                direction = [0 -1 1];
            case 'K'
                direction = [0 0 1];
            case 'L'
                direction = [0 1 1];
            case 'M'
                direction = [0 -1 0];
            case 'N'
                direction = [0 1 0];
            case 'O'
                direction = [0 -1 -1];
            case 'P'
                direction = [0 0 -1];
            case 'Q'
                direction = [0 1 -1];
            case 'R'
                direction = [-1 -1 1];
            case 'S'
                direction = [-1 0 1];
            case 'T'
                direction = [-1 1 1];
            case 'U'
                direction = [-1 -1 0];
            case 'V'
                direction = [-1 0 0];
            case 'W'
                direction = [-1 1 0];
            case 'X'
                direction = [-1 -1 -1];
            case 'Y'
                direction = [-1 0 -1];
            case 'Z'
                direction = [-1 1 -1];
        end
        ends(N,:) = starts(N,:) + direction;
        if N~=length(currentpath)
            starts(N+1,:) = ends(N,:);
        end
    end
    
end

function highlightPath(longestpath,links,Nxyz,dxyz_m,color,spatialFactor)
    dxyz = dxyz_m*spatialFactor.Number;
    clf
    hold on
    for N = 1:1:size(links.start,1)
        plot3([(links.start(N,1)-1)*dxyz(1); (links.end(N,1)-1)*dxyz(1)],[(links.start(N,2)-1)*dxyz(2); (links.end(N,2)-1)*dxyz(2)],[(links.start(N,3)-1)*dxyz(3); (links.end(N,3)-1)*dxyz(3)],'k','LineWidth',1);
    end
    for M = 1:1:size(longestpath.end,1)
        plot3([(longestpath.start(M,1)-1)*dxyz(1); (longestpath.end(M,1)-1)*dxyz(1)],[(longestpath.start(M,2)-1)*dxyz(2); (longestpath.end(M,2)-1)*dxyz(2)],[(longestpath.start(M,3)-1)*dxyz(3); (longestpath.end(M,3)-1)*dxyz(3)],color,'LineWidth',4);
    end
    grid on
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 20;
    xlabel(strcat('$x$-position (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',22,'HorizontalAlignment','left');
    ylabel(strcat('$y$-position (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',22,'HorizontalAlignment','right');
    zlabel(strcat('$z$-position (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',24);
    axis equal
    view(-45,5);
    xlim([0,(Nxyz(1)-1)*dxyz(1)]);
    ylim([0,(Nxyz(2)-1)*dxyz(2)]);    
    zlim([0,(Nxyz(3)-1)*dxyz(3)]);
    set(gcf,'Position',[0,0,1000,1000]);
    set(gcf,'Resize','off');
end