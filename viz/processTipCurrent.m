% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: processTipCurrent.m                                         %
%    Purpose: Identifies the current at the tip of the discharge along    %
%             with the timescale and polarity. Originally embeded         %
%             exclusively within Plot1D_CurrentEstimate.m, but            %
%             applications have since been expanded.                      %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [current, polarity, tau, links, nodes] = processTipCurrent(sims)
    cd ../results/
    EstablishedLinks = load('EstablishedLinks.dat', '-ascii');
    TransportedRho = load('TransportedRhoPos.dat','-ascii');
    TransportedRhoCheck = load('TransportedRhoNeg.dat','-ascii');
    current.polarity = load('TransportedRhoEnd.dat','-ascii');
    cd ../viz/
    
    % Considers the case for a strike depositing charge to the ground/boundary:
    if current.polarity(end)==0 && (EstablishedLinks(end,4)==0 || EstablishedLinks(end,5)==0 || EstablishedLinks(end,6)==0 || (EstablishedLinks(end,4)==(sims.domain.Nx-1) || EstablishedLinks(end,5)==(sims.domain.Ny-1) || EstablishedLinks(end,6)==(sims.domain.Nz-1)))
        current.polarity(end)=-(TransportedRho(end)+TransportedRhoCheck(end));
    % Considers a now-resolved bug that recorded an additional link for ICs:
    elseif abs(current.polarity(end)/current.polarity(end-1))<=(1e-14) && size(current.polarity,1)==(size(EstablishedLinks,1)+1)
        TransportedRho(end)=[];
        TransportedRhoCheck(end)=[];
        current.polarity(end)=[];
    end

    % Convert charge density (nC/m3) into charge (C):
    conversionFactor  = (10^(-9))*sims.domain.dx*sims.domain.dy*sims.domain.dz;
    chargeTransported = conversionFactor.*abs(current.polarity(1:size(EstablishedLinks,1)));
    
    % Isolate polarities:
    polarity.pos = current.polarity;
    polarity.pos(current.polarity>0) = 1;
    polarity.pos(current.polarity<0) = 0;
    
    polarity.neg = current.polarity;
    polarity.neg(current.polarity>0) = 0;
    polarity.neg(current.polarity<0) = 1;
    
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
    
    % Screen print-out of findings:
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

    [links,nodes] = processConnections(EstablishedLinks,tau,polarity,sims);
end