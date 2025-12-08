%% Initialize procedure
close all
clearvars -except sims


cd ../results/

%% Load data
Ez          = load('TotalEfield.dat');
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
phi         = load('TotalPotential.dat');

if isempty(phi) || size(phi,1) == 1
    fprintf('\n*** Plot1D_FieldEvolution.m cannot be executed with current TotalPotential.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_FieldEvolution.m script. ***\n');
end
cd ../viz/

%% Calculate the parameters
z  = ((0:sims.domain.Nz-1)*sims.domain.dz+sims.domain.gnd);      % _m

% Determine magnitude of plotted values:
potentialFactor = checkMagnitude(phi(:));
fieldFactor = checkMagnitude([Ez(:); EthPositive(:); EthNegative(:)]);

NbOfCurves = size(phi);
NbOfCurves = NbOfCurves(1);
CurveStep  = 3;               % Plot every CurveStep
color      = colormap(jet(NbOfCurves));
MaxAbsPhi  = max(max(abs(phi)));
MaxAbsEz   = max(max(abs(Ez))); 

%% Plot figures
figure(1);
set(gcf,'Units','inches','OuterPosition', [20 20 40 20]/6)
% set(gca,'Units','inches','Position', [4 3 15 11]/6)

subplot(121)
hold on
for n=1:NbOfCurves
    if(rem(n,CurveStep)-1==0)
        plot(potentialFactor.Number*phi(n,:),sims.spatialFactor.Number*z,'Color',color(n,:))
    end
end

set(gca,'FontSize',8);
xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);

axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]) ;
box on
grid on

subplot(122)
hold on
for n=1:NbOfCurves
    if(rem(n,CurveStep)-1==0)
        plot(fieldFactor.Number*Ez(n,:),sims.spatialFactor.Number*z,'Color',color(n,:))
    end
end
plot(fieldFactor.Number*EthNegative, sims.spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,sims.spatialFactor.Number*z, 'g--');

set(gca,'FontSize',8);
xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);
axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]);
box on
grid on

%% Multiplots
N = 5;
M = 2;

CurveStep  = floor((NbOfCurves-1)/(M*N-1));
CurveIndex = 1:CurveStep:(NbOfCurves-1);

figure(2);
for n=1:M*N-1
    subplot(M,N,n);
    hold on
    plot(potentialFactor.Number*phi(CurveIndex(n),:),sims.spatialFactor.Number*z,'Color',0*color(CurveIndex(n),:));
    plot(potentialFactor.Number*phi(1,:),sims.spatialFactor.Number*z,'b--')
    hold off
    set(gca,'FontSize',8);
    xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
    ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);
    title(['\phi  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]) ;
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(potentialFactor.Number*phi(NbOfCurves,:),sims.spatialFactor.Number*z,'Color',0*color(NbOfCurves,:));
plot(potentialFactor.Number*phi(1,:),sims.spatialFactor.Number*z,'b--')
set(gca,'FontSize',8);
xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);
title(['\phi  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]) ;
box on
grid on

figure(3);
for n=1:M*N-1
    subplot(M,N,n);
    hold on
    plot(fieldFactor.Number*EthNegative,sims.spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,sims.spatialFactor.Number*z, 'g--');
    plot(fieldFactor.Number*Ez(CurveIndex(n),:),sims.spatialFactor.Number*z,'Color',0*color(CurveIndex(n),:))
    plot(fieldFactor.Number*Ez(1,:),sims.spatialFactor.Number*z,'r--')
    hold off
    set(gca,'FontSize',8);
    xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
    ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);
    title(['E_z  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]);
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(fieldFactor.Number*EthNegative,sims.spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,sims.spatialFactor.Number*z, 'g--');
plot(fieldFactor.Number*Ez(NbOfCurves,:),sims.spatialFactor.Number*z,'Color',0*color(NbOfCurves,:))
plot(fieldFactor.Number*Ez(1,:),sims.spatialFactor.Number*z,'r--')
hold off
set(gca,'FontSize',8);
xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',10);
title(['E_z  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)]);
box on
grid on
clear n m M N CurveIndex CurveStep NbOfCurves