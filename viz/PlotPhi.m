close all
clear all
beep  off
clc

d     = load('d.dat');
N     = load('N.dat');
z_gnd = load('z_gnd.dat');
nn = 1;

for ii=211650:50:211700
    Ez_tmp    = load(['Ez2d',num2str(ii),'.dat']);
    Ez(nn,:)  = Ez_tmp(1,:);
    nn = nn+1;
end
    
EthPositive = load('Ec.dat')*2.10/2.16;
EthNegative = -load('Ec.dat')*2.10/2.16;

%%
close all
clc
z = z_gnd + (0:N(2)-1)*d(2);
plot(Ez(:,:),z,'-',EthPositive,z,'g-.',EthNegative,z,'g-.');


stop

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
EzNum       = load('EnumBF.dat');
phiNum      = load('phiNumBF.dat');
z_gnd       = load('z_gnd.dat');
E2D         = load('E2D.dat');
phi2D       = load('phi2D.dat');
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
Ez          = load('Ez2d120900.dat');

dy = dxyz(2);
dz = dxyz(3);
Ny = Nxyz(2);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
Y = (0:dy:(Ny-1)*dy);
Z = (0:dz:(Nz-1)*dz);

figure;
hold on
plot(...
    EzNum,Z+z_gnd,'b',...
    1/1000*phiNum,Z+z_gnd,'r',...
    Ez(1,1:length(Z)),Z+z_gnd,'k',...
    EthNegative,Z+z_gnd,'g',...
    EthPositive,Z+z_gnd,'g'...
    );
plot([-max(max(abs(EzNum)),max(1/1000*abs(phiNum))) max(max(abs(EzNum)),max(1/1000*abs(phiNum)))],[z_gnd,z_gnd],'k');
hold off
% axis([1e+2 1e+7 0 40e+3]);
axis([...
    -max(max(abs(EzNum)),max(1/1000*abs(phiNum))) ...
    max(max(abs(EzNum)),max(1/1000*abs(phiNum))) ...
    0 Lz+z_gnd]);
legend('E_z','\phi','E_z from 2-D');
xlabel('E_z (kV/cm) & \phi (V)');
ylabel('Altitude (m)');
grid on

% figure;
% imagesc(Z+z_gnd,Y,E2D);
% ylabel('y-axis (m)');
% xlabel('z-axis (m)');
% title('|E| in the plane (i_c , : , :)');
stop
figure;
imagesc(Z+z_gnd,Y,phi2D);
ylabel('y-axis (m)');
xlabel('z-axis (m)');
title('\phi  in the plane (i_c , : , :)');
colorbar;