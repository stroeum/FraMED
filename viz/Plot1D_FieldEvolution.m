%% Initialize procedure
close all
clearvars -except sims


cd ../results/

%% Load data
dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
Ez          = load('TotalEfield.dat');
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
phi         = load('TotalPotential.dat');
z_gnd       = load('z_gnd.dat');

if isempty(phi) || size(phi,1) == 1
    fprintf('\n*** Plot1D_FieldEvolution.m cannot be executed with current TotalPotential.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_FieldEvolution.m script. ***\n');
end
cd ../viz/

%% Calculate the parameters
dz = dxyz(3);                  % _m
Nz = Nxyz(3);                  % _
Lz = (Nz-1)*dz;                % _m
z  = (z_gnd+(0:Nz-1)*dz);      % _m
clear dxyz Nxyz

% Determine magnitude of plotted values:
potentialFactor = checkMagnitude(phi(:));
spatialFactor = checkMagnitude(z(:));
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
        plot(potentialFactor.Number*phi(n,:),spatialFactor.Number*z,'Color',color(n,:))
    end
end

set(gca,'FontSize',8);
xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);

axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi spatialFactor.Number*min(z) spatialFactor.Number*max(z)]) ;
box on
grid on

subplot(122)
hold on
for n=1:NbOfCurves
    if(rem(n,CurveStep)-1==0)
        plot(fieldFactor.Number*Ez(n,:),spatialFactor.Number*z,'Color',color(n,:))
    end
end
plot(fieldFactor.Number*EthNegative, spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,spatialFactor.Number*z, 'g--');

set(gca,'FontSize',8);
xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);
axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz spatialFactor.Number*min(z) spatialFactor.Number*max(z)]);
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
    plot(potentialFactor.Number*phi(CurveIndex(n),:),spatialFactor.Number*z,'Color',0*color(CurveIndex(n),:));
    plot(potentialFactor.Number*phi(1,:),spatialFactor.Number*z,'b--')
    hold off
    set(gca,'FontSize',8);
    xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
    ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);
    title(['\phi  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi spatialFactor.Number*min(z) spatialFactor.Number*max(z)]) ;
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(potentialFactor.Number*phi(NbOfCurves,:),spatialFactor.Number*z,'Color',0*color(NbOfCurves,:));
plot(potentialFactor.Number*phi(1,:),spatialFactor.Number*z,'b--')
set(gca,'FontSize',8);
xlabel(strcat('\phi (',potentialFactor.Unit,'V)'),'FontSize',10);
ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);
title(['\phi  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-potentialFactor.Number*MaxAbsPhi potentialFactor.Number*MaxAbsPhi spatialFactor.Number*min(z) spatialFactor.Number*max(z)]) ;
box on
grid on

figure(3);
for n=1:M*N-1
    subplot(M,N,n);
    hold on
    plot(fieldFactor.Number*EthNegative,spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,spatialFactor.Number*z, 'g--');
    plot(fieldFactor.Number*Ez(CurveIndex(n),:),spatialFactor.Number*z,'Color',0*color(CurveIndex(n),:))
    plot(fieldFactor.Number*Ez(1,:),spatialFactor.Number*z,'r--')
    hold off
    set(gca,'FontSize',8);
    xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
    ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);
    title(['E_z  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz spatialFactor.Number*min(z) spatialFactor.Number*max(z)]);
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(fieldFactor.Number*EthNegative,spatialFactor.Number*z, 'g--',fieldFactor.Number*EthPositive,spatialFactor.Number*z, 'g--');
plot(fieldFactor.Number*Ez(NbOfCurves,:),spatialFactor.Number*z,'Color',0*color(NbOfCurves,:))
plot(fieldFactor.Number*Ez(1,:),spatialFactor.Number*z,'r--')
hold off
set(gca,'FontSize',8);
xlabel(strcat('E_z (',fieldFactor.Unit,'V/m)'),'FontSize',10);
ylabel(strcat('z (',spatialFactor.Unit,'m)'),'FontSize',10);
title(['E_z  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-fieldFactor.Number*MaxAbsEz fieldFactor.Number*MaxAbsEz spatialFactor.Number*min(z) spatialFactor.Number*max(z)]);
box on
grid on
clear n m M N CurveIndex CurveStep NbOfCurves