%% Initialization
close all
clearvars
clc

cd ../results

%% Load Data
dxyz   = load('dxyz.dat');
Nxyz   = load('Nxyz.dat');
z_gnd  = load('z_gnd.dat')*1e-3;
phi.XZ = load('phi2D.dat');
phi.YZ = load('phi2D.dat');

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
contour(X,Z,phi.XZ,60);
contourf(X,Z,phi.XZ,60,'LineColor','none');
caxis([-max(max(abs(phi.XZ))) max(max(abs(phi.XZ)))])
% imagesc(Y,Z,phi.XZ');
axis image
axis xy
ylabel('z-axis (m)','FontSize',12);
xlabel('x-axis (m)','FontSize',12);
set(gca,'XMinorTick','on','YMinorTick','on');
title('\phi (V)','FontSize',12);
colorbar;

subplot(122);
contour(Y,Z,phi.YZ,50);
contourf(Y,Z,phi.YZ,60,'LineColor','none');
caxis([-max(max(abs(phi.YZ))) max(max(abs(phi.YZ)))])
% imagesc(Y,Z,phi.YZ');
axis image
axis xy
ylabel('z-axis (m)','FontSize',12);
xlabel('y-axis (m)','FontSize',12);
set(gca,'XMinorTick','on','YMinorTick','on');
title('\phi (V)','FontSize',12);
colorbar;

cd ../viz