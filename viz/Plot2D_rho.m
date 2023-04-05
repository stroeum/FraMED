%% Initialization
close all
clearvars -except sims


cd ../results

%% Load Data
dxyz   = load('dxyz.dat');
Nxyz   = load('Nxyz.dat');
z_gnd  = load('z_gnd.dat')*1e-3;
rho.XZ = load('rhoAmbXZ.dat');
rho.YZ = load('rhoAmbYZ.dat');

d.x    = dxyz(1)*1e-3;
d.y    = dxyz(2)*1e-3;
d.z    = dxyz(3)*1e-3;
N.x    = Nxyz(1);
N.y    = Nxyz(2);
N.z    = Nxyz(3);
clear dxyz Nxyz

X = (0:d.x:(N.x-1)*d.x);
Y = (0:d.y:(N.y-1)*d.y);
Z = (0:d.z:(N.z-1)*d.z)+z_gnd;

%% Plot
figure;
subplot(121);
contour(X,Z,rho.XZ,60);
contourf(X,Z,rho.XZ,60,'LineColor','none');
caxis([-max(max(abs(rho.XZ))) max(max(abs(rho.XZ)))])
% imagesc(Y,Z,rho.XZ');
axis image
axis xy
ylabel('z-axis (m)','FontSize',12);
xlabel('x-axis (m)','FontSize',12);
set(gca,'XMinorTick','on','YMinorTick','on');
title('\rho_{amb} (nC/m^3)','FontSize',12);
colorbar;

subplot(122);
contour(Y,Z,rho.YZ,50);
contourf(Y,Z,rho.YZ,60,'LineColor','none');
caxis([-max(max(abs(rho.YZ))) max(max(abs(rho.YZ)))])
% imagesc(Y,Z,rho.YZ');
axis image
axis xy
ylabel('z-axis (m)','FontSize',12);
xlabel('y-axis (m)','FontSize',12);
set(gca,'XMinorTick','on','YMinorTick','on');
title('\rho_{amb} (nC/m^3)','FontSize',12);
colorbar;

cd ../viz