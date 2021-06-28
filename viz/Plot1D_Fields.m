function []=Plot1D_Fields()
clear all
close all
clc

%% Compile and Run C++ code
% cd ..
% !make
% !./main
% cd results

%% Plot results

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
phiNumBF    = load('phiNumBF.dat');
%phiNumAF    = load('../results/phiNumAF.dat');
phiAn       = load('phiAnalytical.dat');
%phiAn       = load('../results copy/phiNumAF.dat');
z_gnd       = load('z_gnd.dat');

dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz

linewidth = 1;
figure;
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(phiNumBF(:)*1e-6, z*1e-3, 'b-','LineWidth',linewidth)
%plot(phiNumAF(:)*1e-6, z*1e-3, 'b--','LineWidth',linewidth)
plot(phiAn(:)*1e-6, z*1e-3, 'r--','LineWidth',linewidth)
plot([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6)],[z_gnd,z_gnd]*1e-3,'k','LineWidth',linewidth);
hold off
axis([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6) 0 (Lz+z_gnd)*1e-3]);
xlabel('\phi (MV)','FontSize',16);
ylabel('Altitude (km)','FontSize',16);
set(gca,'FontSize',14);
set(gca,'XMinorTick','on','YMinorTick','on')
box on
grid on

end