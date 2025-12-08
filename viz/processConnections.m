% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: processConnections.m                                        %
%    Purpose: Identifies the link structures connections and begins to    %
%             process the current and charge densities. Originally        %
%             embeded exclusively within Plot1D_CurrentEstimate.m, but    %
%             applications have since been expanded.                      %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [links,nodes] = processConnections(EstablishedLinks,tau,polarity,sims)
    % Convert from nC/m3 to C:
    conversionFactor  = (10^(-9))*sims.domain.dx*sims.domain.dy*sims.domain.dz;
    
    links.start    = EstablishedLinks(:,1:3);
    links.end      = EstablishedLinks(:,4:6);
    links.distance = EstablishedLinks(:,7);
    links.initiate = links.start(1,:);
    links.steps    = 1:1:size(EstablishedLinks,1);
    links.height   = EstablishedLinks(:,6)*sims.domain.dz;
    
    % Assign important values regarding the links:
    sims.spatialFactor  = checkMagnitude(links.height);

    % Identifies total branch lengths to predict lower limit of timescale:
    temp = load("../results/NodeDeltaRho.dat",'-ascii');
    temp2 = load("../results/NodeRho.dat",'-ascii');
    
    % Determines various path branches and distances along each branch:
    statusBar = waitbar(0,strcat("(",num2str(0,'%.0f'),"%) Processing ",num2str(0)," out of ",num2str(size(EstablishedLinks,1))," links."));
    links.total = size(EstablishedLinks,1);
    
    links.connections.num = zeros([links.total,(links.total+1)]);
    links.connections.nodes = zeros([links.total,1]);
    links.connections.deltas = zeros([links.total,(links.total)+1]);
    links.connections.currents = zeros([links.total,(links.total)+1]);
    nodes.rho.deltas = NaN([links.total,(links.total+1)]);
    nodes.rho.values = NaN([links.total,(links.total+1)]);
    
    tracker = 1;
    for N = 1:1:links.total
        waitbar((N/links.total),statusBar,strcat("(",num2str(100*(N-1)/links.total,'%.0f'),"%) Processing ",num2str(N)," out of ",num2str(size(EstablishedLinks,1))," links."));
        if N == 1
            links.startnode = 1;
            links.endnode = 2;
            links.connections.num(N,1) = 1;
            links.connections.num(N,2) = 1;
        else
            links.endnode = N+1;
            links.connections.num(N,:) = links.connections.num((N-1),:);
            links.connections.num(N,N+1) = 1;
        end
        
        % If the starting node is the initiation point:
        if links.start(N,:) == links.initiate
            % Classify path for this particular link:
            links.path(N) = string(characterizePropagation(links,N));
            % Total distance along this path/branch:
            links.pathdistance(N) = links.distance(N);
            % Track total distance of all links:
            if N == 1
                links.totaldistance(N) = links.distance(N);
            else
                links.totaldistance(N) = links.totaldistance(N-1)+links.distance(N);
                links.connections.num(N,1) = links.connections.num(N,1) + 1;
            end
            links.connections.nodes(N) = 0;
            links.mintime(N) = tau.case(N);
            links.maxtime(N) = tau.case(N);
        % Otherwise, if the start node is not at the initiation point...
        else
            % Determine which path/branch this link is an extension of:
            match = find(ismember(links.end(1:(N-1),:),links.start(N,:),'row'));
            links.connections.num(N,match+1) = links.connections.num(N,match+1) + 1;
            links.connections.nodes(N) = match;
            % Classify path for this particular link:
            links.path(N) = string(strcat(string(links.path(match)),string(characterizePropagation(links,N))));
            reshape(links.path,[N,1]);
            % Total distance along this path/branch:
            links.pathdistance(N) = links.pathdistance(match) + links.distance(N);
            % Track total distance of all links:
            links.totaldistance(N) = links.totaldistance(N-1)+links.distance(N);
            links.mintime(N) = links.mintime(match) + tau.case(N);
            links.maxtime(N) = links.maxtime(N-1) + tau.case(N);
        end
        
        % Convert changes in charge density (nC/m3) into current (C/s):
        nodes.rho.deltas(N,1:N+1) = temp(tracker:tracker+N).*conversionFactor./tau.case(N);
    
        % Convert local charge density array to matrix (nC/m3):
        nodes.rho.values(N,1:N+1) = temp2(tracker:tracker+N);
        tracker = tracker+N+1;
    end
    nodes.rho.intensity = abs(nodes.rho.values)./max(max(abs(nodes.rho.values)));
    close(statusBar);
    [links.maxdistance, links.maxindex] = max(links.pathdistance);
    links.timescale.max = sum(tau.case);
    links.timescale.min = links.maxdistance/(polarity.neg(links.maxindex)*sims.vprop.neg + polarity.pos(links.maxindex)*sims.vprop.pos);
end