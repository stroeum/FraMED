%% Initialize procedure
close all
clearvars -except sims


cd ../results/

%% Load data
dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
Ez          = load('TotalEfield.dat')*1e-5;
EthPositive = load('EthPositive.dat')*1e-5;
EthNegative = load('EthNegative.dat')*1e-5;
phi         = load('TotalPotential.dat')*1e-6;
z_gnd       = load('z_gnd.dat');

if isempty(phi)
    fprintf('\n*** Plot1D_FieldEvolution.m cannot be executed with current TotalPotential.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_FieldEvolution.m script. ***\n');
end

%% Calculate the parameters
dz = dxyz(3);                  % _m
Nz = Nxyz(3);                  % _
Lz = (Nz-1)*dz*1e-3;           % _km
z  = (z_gnd+(0:Nz-1)*dz)*1e-3; % _km
clear dxyz Nxyz

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
        plot(phi(n,:),z,'Color',color(n,:))
    end
end
set(gca,'FontSize',8);
xlabel('\phi (MV)','FontSize',10);
ylabel('z(km)','FontSize',10);

axis([-MaxAbsPhi MaxAbsPhi min(z) max(z)]) ;
box on
grid on

subplot(122)
hold on
for n=1:NbOfCurves
    if(rem(n,CurveStep)-1==0)
        plot(Ez(n,:),z,'Color',color(n,:))
    end
end
plot(EthNegative,z, 'g--',EthPositive,z, 'g--');

set(gca,'FontSize',8);
xlabel('E_z (kV/cm)','FontSize',10);
ylabel('z(km)','FontSize',10);
axis([-MaxAbsEz MaxAbsEz min(z) max(z)]);
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
    plot(phi(CurveIndex(n),:),z,'Color',0*color(CurveIndex(n),:));
    plot(phi(1,:),z,'b--')
    hold off
    set(gca,'FontSize',8);
    xlabel('\phi (MV)','FontSize',10);
    ylabel('z(km)','FontSize',10);
    title(['\phi  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-MaxAbsPhi MaxAbsPhi min(z) max(z)]) ;
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(phi(NbOfCurves,:),z,'Color',0*color(NbOfCurves,:));
plot(phi(1,:),z,'b--')
set(gca,'FontSize',8);
xlabel('\phi (MV)','FontSize',10);
ylabel('z(km)','FontSize',10);
title(['\phi  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-MaxAbsPhi MaxAbsPhi min(z) max(z)]);
box on
grid on

figure(3);
for n=1:M*N-1
    subplot(M,N,n);
    hold on
    plot(EthNegative,z, 'g--',EthPositive,z, 'g--');
    plot(Ez(CurveIndex(n),:),z,'Color',0*color(CurveIndex(n),:))
    plot(Ez(1,:),z,'r--')
    hold off
    set(gca,'FontSize',8);
    xlabel('E_z (kV/cm)','FontSize',10);
    ylabel('z(km)','FontSize',10);
    title(['E_z  after ',int2str(CurveIndex(n)-1),' steps'],'FontSize',5);
    axis([-MaxAbsEz MaxAbsEz min(z) max(z)]);
    box on
    grid on
end
subplot(M,N,M*N);
hold on
plot(EthNegative,z, 'g--',EthPositive,z, 'g--');
plot(Ez(NbOfCurves,:),z,'Color',0*color(NbOfCurves,:))
plot(Ez(1,:),z,'r--')
hold off
set(gca,'FontSize',8);
xlabel('E_z (kV/cm)','FontSize',10);
ylabel('z(km)','FontSize',10);
title(['E_z  after ',int2str(NbOfCurves-1),' steps'],'FontSize',5);
axis([-MaxAbsEz MaxAbsEz min(z) max(z)]);
box on
grid on
clear n m M N CurveIndex CurveStep NbOfCurves

cd ../viz/