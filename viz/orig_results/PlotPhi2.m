% load('EzNum.dat')
% load('EzAna.dat')
% load('dxyz.dat');
% z=[0:size(EzNum)-1]*dxyz(3);
% plot(EzNum,z, 'r', EzAna,z,'b')
% semilogx(EzNum,z, 'r', EzAna,z,'b')

% load('EnumBF.dat')
% load('dxyz.dat');
% z=[0:size(EnumBF)-1]*dxyz(3);
% plot(EnumBF,z, 'r')

close all
clear all
clc

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
z_gnd       = load('z_gnd.dat');
phiCa       = load('phiC1Da.dat');
phiCb       = load('phiC1Db.dat');
phiA        = load('phiA1D.dat');
phiC2D      = load('phiC2D.dat');
phiA2D      = load('phiA2D.dat');

dy = dxyz(2);
dz = dxyz(3);
Ny = Nxyz(2);
Nz = Nxyz(3);
Y = (0:dy:(Ny-1)*dy);
Z = (0:dz:(Nz-1)*dz);

figure;
subplot(211)
plot(...
    phiCa       ,Z+z_gnd,'b',...
    phiA        ,Z+z_gnd,'r',...
    phiA +phiCa ,Z+z_gnd,'g'...
    );
xlabel('\phi (V)');
ylabel('Altitude (m)');
legend('\phi_{channel}','\phi_{ambient}','\phi_{total}');
title('\phi  before SOR adjustment');
grid on

subplot(212)
plot(...
    phiCb       ,Z+z_gnd,'b',...
    phiA        ,Z+z_gnd,'r',...
    phiA +phiCb ,Z+z_gnd,'g'...
);
xlabel('\phi (V)');
ylabel('Altitude (m)');
legend('\phi_{channel}','\phi_{ambient}','\phi_{total}');
title('\phi  after SOR adjustment');
grid on

figure;
subplot(211)
imagesc(Z+z_gnd,Y,phiC2D);
ylabel('y-axis (m)');
xlabel('z-axis (m)');
title('\phi_{channel}');
colorbar;
subplot(212)
imagesc(Z+z_gnd,Y,phiA2D);
ylabel('y-axis (m)');
xlabel('z-axis (m)');
title('\phi_{channel}');
colorbar;

figure;
imagesc(Z+z_gnd,Y,phiC2D+phiA2D);
ylabel('y-axis (m)');
xlabel('z-axis (m)');
title('\phi_{total}');
colorbar;