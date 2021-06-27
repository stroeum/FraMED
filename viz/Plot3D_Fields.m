function []=Plot3D_Fields()
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
EzNum       = load('EzNum.dat');
% phiNum      = load('phiNum.dat');
% Einitiation = load('Einitiation.dat');
% EthPositive = load('EthPositive.dat');
% EthNegative = load('EthNegative.dat');
z_gnd       = load('z_gnd.dat');

gnd_color = [.5 .5 .5];

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

clear Nxyz
clear dxyz
z  = (0:Nz-1)*dz; 

% plot(...
%     EzNum(:), z,'r-.',...
%     1/1000*phiNum(:), z, 'b-'...
%     );
% 
% legend('E_z','\phi');
% xlabel('E_z (V/m) & \phi (V)');
% ylabel('Altitude (m)');
% grid on

x             = EzNum/max(max(max(abs(EzNum))))*1e3;
y             = z;
[theta,rho,z] = cart2pol(x, zeros(size(x,1), size(x,2)), y); 
r             = 10;
theta         = repmat( linspace(0,2*pi,r)', 1, length(rho(:)) );
rho           = repmat(rho(:).', r, 1);
z             = repmat(z(:).', r, 1);
[X,Y,Z]       = pol2cart(theta,rho,z);
% waterfall(X,Y,Z);
hold on
mesh(X+Lx/2,Y+Ly/2,Z+z_gnd,'FaceColor','none');
Xp = [Lx*3/4 Lx/4 Lx/4 Lx*3/4]; 
Yp = [Ly*3/4 Ly*3/4 Ly/4 Ly/4]; 
Zp = [z_gnd z_gnd z_gnd z_gnd];
patch(Xp, Yp, Zp, z_gnd,'FaceColor',gnd_color);
hold off
axis([Lx/4 Lx*3/4 Ly/4 Ly*3/4 0 Lz+z_gnd])
colormap([0 0 0])
end