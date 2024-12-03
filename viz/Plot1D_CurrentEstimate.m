close all
clearvars -except sims 

fprintf('\n*** Executing Plot1D_CurrentEstimate.m script. ***\n');  
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    sims = specifySimDetails();
end

%% Load data from results folder:
cd ../results/
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
load('EstablishedLinks.dat', '-ascii');
load('TransportedRho.dat','-ascii');
cd ../viz/

%% Calculate current:
% Convert charge density (nC/m3) into charge (C):
conversionFactor  = (10^(-9))*dxyz(1)*dxyz(2)*dxyz(3);
chargeTransported = conversionFactor.*TransportedRho;

% Declared constants:
vprop.leader   = 4.4*(10^5); % propagation speed for leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
vprop.streamer = 3.0*(10^7); % propagation speed for streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9

% Calculates instantaneous timescale on a per-link basis:
tau.leader.avg   = EstablishedLinks(:,7)./vprop.leader;
tau.streamer.avg = EstablishedLinks(:,7)./vprop.streamer;

% Calculates the current (units of kiloAmps):
current.leader.avg   = 0.001*chargeTransported./tau.leader.avg;
current.streamer.avg = 0.001*chargeTransported./tau.streamer.avg;

% Assigns the correct values for the system:
if strcmp(sims.objectType,'Streamer')
    vprop.case   = vprop.streamer;
    tau.case     = tau.streamer.avg;
    current.case = current.streamer.avg;
elseif strcmp(sims.objectType,'Leader')
    vprop.case   = vprop.leader;
    tau.case     = tau.leader.avg;
    current.case = current.leader.avg;
end

%% Identifies total branch lengths to predict lower limit of timescale:
% Assign important values:
links.start    = EstablishedLinks(:,1:3);
links.end      = EstablishedLinks(:,4:6);
links.distance = EstablishedLinks(:,7);
links.initiate = links.start(1,:);
links.steps    = 1:1:size(EstablishedLinks,1);

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
        % Track instanteneous current:
        links.current(N) = vprop.case*(links.pathcharge(N))/links.pathdistance(N);

    end
end
close(statusBar);
[links.maxdistance, maxindex] = max(links.pathdistance);
links.timescale.max = sum(tau.case);
links.timescale.min = links.maxdistance/vprop.case;

% Highlights longest path (or path to ground, if applicable):
[longestpath.start, longestpath.end] = string2link(links.path(maxindex),links.initiate);
highlightPath(longestpath,links,Nxyz,dxyz,'r')
exportgraphics(gcf,[sims.pathPNGs,'/Path-Longest_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);  % Before discharge
if links.end(end,3)==0
    [groundpath.start, groundpath.end] = string2link(links.path(end),links.initiate);
    highlightPath(groundpath,links,Nxyz,dxyz,'b')
    exportgraphics(gcf,[sims.pathPNGs,'/Path-Ground_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);  % Before discharge
end

%% Screen print-out of findings:
[current.min.val, current.min.index] = min(current.case);
[current.max.val, current.max.index] = max(current.case);
current.avg.val = mean(current.case);
current.avg.std = std(current.case);
current.partial.val = mean(current.case(1:(end-1)));
current.partial.std = std(current.case(1:(end-1)));
fprintf(strcat("\n\t",sims.objectType," propagation speed approximated as ",num2str(vprop.case/1000,'%.0f')," km/s."))
fprintf(strcat("\n\tDischarge timescale estimated to be between ",num2str(1000*links.timescale.min,'%.1f')," - ",num2str(1000*links.timescale.max,'%.1f')," ms."))
fprintf(strcat("\n\tLongest branch traveled ",num2str(links.maxdistance*0.001,'%.1f')," km from the initiation point."))
if strcmp(sims.objectType,'Leader')
    fprintf(strcat("\n\tAverage currents (based on charge transfer):"));
    fprintf(strcat("\n\t\t",num2str(current.avg.val,5),"\tkA\t(mean, including all links,  sigma = ",num2str(current.avg.std,5)," kA)"));
    fprintf(strcat("\n\t\t",num2str(current.partial.val,5),"\tkA\t(mean, excluding final link, sigma = ",num2str(current.partial.std,5)," kA)"));
    fprintf(strcat("\n\t\t",num2str(current.max.val,5),"\tkA\t(maximum, link ",num2str(links.steps(current.max.index)),")"));
    fprintf(strcat("\n\t\t",num2str(current.min.val,5),"\tkA\t(minimum, link ",num2str(links.steps(current.min.index)),")\n"));
else
    fprintf(strcat("\n\tAverage currents (based on charge transfer):"));
    fprintf(strcat("\n\t\t",num2str(1000*current.avg.val,5),"\tA\t(mean, including all links,  sigma = ",num2str(1000*current.avg.std,5)," A)"));
    fprintf(strcat("\n\t\t",num2str(1000*current.partial.val,5),"\tA\t(mean, excluding final link, sigma = ",num2str(1000*current.partial.std,5)," A)"));
    fprintf(strcat("\n\t\t",num2str(1000*current.max.val,5),"\tA\t(maximum, link ",num2str(links.steps(current.max.index)),")"));
    fprintf(strcat("\n\t\t",num2str(1000*current.min.val,5),"\tA\t(minimum, link ",num2str(links.steps(current.min.index)),")\n"));
end
if links.end(end,3)==0
    fprintf(strcat("\t\t",num2str(0.001*chargeTransported(end)./(links.pathdistance(end)/vprop.case),5),"\tkA\t(initiation-to-ground path)\n"))
end

%% Plot of instantaneous current values:
clf
hold on; box on;
scatter(links.steps,current.case,25,'r.','DisplayName','Instantaneous Current: $I_l = \frac{\Delta  Q_l}{\Delta t_l}$')
yline(current.max.val,'k:','DisplayName',strcat("Maximum Current: ",num2str(current.max.val,3)," kA (link ",num2str(links.steps(current.max.index)),")"),'LineWidth',2);
if links.end(end,3)==0
    yline(current.avg.val,'k--','DisplayName',strcat("Average Current: ",num2str(current.partial.val,3)," $\pm$ ",num2str(current.partial.std,3)," kA"),'LineWidth',2);
else
    yline(current.avg.val,'k--','DisplayName',strcat("Average Current: ",num2str(current.avg.val,3)," $\pm$ ",num2str(current.avg.std,3)," kA"),'LineWidth',2);
end
yline(current.min.val,'k-','DisplayName',strcat("Minimum Current: ",num2str(current.min.val,3)," kA (link ",num2str(links.steps(current.min.index)),")"),'LineWidth',2);
xlim([0,links.steps(end)]);
legend('Interpreter','latex','FontSize',16,'location','best','box','off')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
title(strcat("Current Estimate for ",sims.objectType," Discharge on ",sims.objectName),'Interpreter','latex','FontSize',28);
xlabel('Propagation Link $l$','FontSize',20,'Interpreter','latex');
ylabel('Current $I$ (kA)','FontSize',20,'Interpreter','latex');
grid off; hold off;
set(gcf,'Position',[0,0,1200,800]);
exportgraphics(gcf,[sims.pathPNGs,'/CurrentEstimate_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);  

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

function highlightPath(longestpath,links,Nxyz,dxyz_m,color)
    dxyz = dxyz_m/1000;
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
    xlabel('$x$-position (km)','Interpreter','latex','FontSize',22,'HorizontalAlignment','left');
    ylabel('$y$-position (km)','Interpreter','latex','FontSize',22,'HorizontalAlignment','right');
    zlabel('$z$-position (km)','Interpreter','latex','FontSize',24);
    axis equal
    view(-45,9);
    xlim([0,(Nxyz(1)-1)*dxyz(1)]);
    ylim([0,(Nxyz(2)-1)*dxyz(2)]);    
    zlim([0,(Nxyz(3)-1)*dxyz(3)]);
    set(gcf,'Position',[0,0,1000,1000]);
    set(gcf,'Resize','off');
end