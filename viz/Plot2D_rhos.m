close all
clearvars
clc

cd ../results

%% Load datas

dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
rhos2D       = load('rhoAmbYZ0.dat')*1e-9;
z_gnd        = load('z_gnd.dat');
ChargeLayers = load('ChargeLayers.dat');

dx = dxyz(1);               % _m
dy = dxyz(2);               % _m
dz = dxyz(3);               % _m

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

Lx = (Nx-1)*dx*1e-3;         % _m
Ly = (Ny-1)*dy*1e-3;         % _m
Lz = (Nz-1)*dz*1e-3;         % _m

clear dxyz Nxyz

fprintf('SOR_data loaded and saved\n')

%% Derive data for plotting
x        = (0:Nx-1)'*dx*1e-3;
y        = (0:Ny-1)'*dy*1e-3;
z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
% [x,y]    = meshgrid(x,y);

%% Plot
figure;
hold on
set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/3)
axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))]);

imagesc(x,y,rhos2D);
colorbar;

hold off
xlabel('x (km)','FontSize',12);
ylabel('y (km)','FontSize',12);
%title('Projection of E-field lines','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',10);
axis square
box on

cd ../viz/