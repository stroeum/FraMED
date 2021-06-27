close all
clear all
clc

cd ..
cd resultsClosedBC/
dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
Ex2D         = load('Ex2D.dat');
Ey2D         = load('Ey2D.dat');
Ez2D         = load('Ez2D.dat');
phi2D        = load('phi2D.dat');
z_gnd        = load('z_gnd.dat');

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

y  = (0:Ny-1)*dy*1e-3;
z  = ((0:Nz-1)*dz + z_gnd)*1e-3;
ko = find( z == 9.8);
Eo = ((Ex2D.^2+Ey2D.^2+Ez2D.^2)).^.5*1e-5;

y_CBC  = y;
Eo_CBC = Eo(:,ko);

cd ..
cd resultsOpenBC/
dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
Ex2D         = load('Ex2D.dat');
Ey2D         = load('Ey2D.dat');
Ez2D         = load('Ez2D.dat');
phi2D        = load('phi2D.dat');
z_gnd        = load('z_gnd.dat');
cd ..

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

y  = (0:Ny-1)*dy*1e-3;
z  = ((0:Nz-1)*dz + z_gnd)*1e-3;
ko = find( z == 9.8);
Eo = ((Ex2D.^2+Ey2D.^2+Ez2D.^2)).^.5*1e-5;

y_OBC  = y;
Eo_OBC = Eo(:,ko);

cd ..

minE = .5*min(min(Eo_OBC),min(Eo_CBC));
maxE = 1.5*max(max(Eo_OBC),max(Eo_CBC));
miny = min(y_OBC);
maxy = max(y_OBC);

semilogy(y_CBC,Eo_CBC,y_OBC,Eo_OBC);
axis([miny maxy minE maxE]);
legend('Closed BC','Open BC','Location','South');
xlabel('y (km)','FontSize',12);
ylabel('|E| (kV/cm)','FontSize',12);
title('E-field magnitude (kV/cm) @ Alt. 9.8 km','FontSize',12,'FontWeight','bold');
grid on