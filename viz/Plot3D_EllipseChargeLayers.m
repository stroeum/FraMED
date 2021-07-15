close all
clearvars
clc

cd ../results/
dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
z_gnd        = load('z_gnd.dat');
ChargeLayers = load('ChargeLayers.dat');
cd ../viz/

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

clear Nxyz dxyz

NbVer      = 25;
NbRad      = 25;
NbOrthoRad = 25;
Color      = ['r' 'r' 'b' 'b'];

hold on
for n = 1:4
    Ellipse(ChargeLayers(n,2)*1e-3,ChargeLayers(n,3)*1e-3,ChargeLayers(n,4)*1e-3,...
        ChargeLayers(n,5)*1e-3,ChargeLayers(n,6)*1e-3,ChargeLayers(n,7)*1e-3,...
        Color(n),NbVer,NbRad,NbOrthoRad);
end
set(gca,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[1 1 1]);
axis([0 Lx 0 Ly z_gnd z_gnd+Lz]*1e-3);
xlabel('x-axis (km)');
ylabel('y-axis (km)');
zlabel('z-axis (km)');
axis on
box on
grid on
view([1 0 0])
hold off
