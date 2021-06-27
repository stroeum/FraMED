function []=Plot1D_Fields()
clear all
close all
clc

%% Compile and Run C++ code
% cd ..
% !rm -rf results/*.dat results/*.avi results/*.mat
% !make
% !./main
% cd results

%% Plot results

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
EzNumBF     = load('EnumBF.dat');
phiNumBF    = load('phiNumBF.dat');
EzNumAF     = load('EnumAF.dat');
phiNumAF    = load('phiNumAF.dat');
% Einitiation = load('Einitiation.dat');
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
z_gnd       = load('z_gnd.dat');

dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz

figure;
subplot(121)
hold on
plot(...
    EzNumBF(:),z*1e-3,'r-.',...
    1/1000*phiNumBF(:), z*1e-3, 'b-',...
    EthNegative(:),z*1e-3, 'g--',...
    EthPositive(:),z*1e-3, 'g--'...
    );
plot([min(EthNegative) max(EthPositive)],[z_gnd,z_gnd]*1e-3,'k');
% semilogx(...
%     abs(EzNum(:)), [0:size(EzNum)-1]*dz,'r-',...
%     abs(Einitiation(:)),[0:size(Einitiation)-1]*dz, 'b--',...
%     abs(EthNegative(:)),[0:size(EthNegative)-1]*dz, 'g--',...
%     abs(EthPositive(:)),[0:size(EthPositive)-1]*dz, 'g--'...
%     );

axis([min(EthNegative) max(EthPositive) 0 (Lz+z_gnd)*1e-3]);
legend('E_z [BF]','\phi [BF]','E_{th}^\pm');
xlabel('E_z (V/m) & \phi (kV)');
ylabel('Altitude (km)');
title('E_z and \phi  before the flash [BF]');
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
box on
hold off

subplot(122)
hold on
plot(...
    EzNumAF(:),z*1e-3,'r-.',...
    1/1000*phiNumAF(:), z*1e-3, 'b-',...
    EthNegative(:),z*1e-3, 'g--',...
    EthPositive(:),z*1e-3, 'g--'...
    );
plot([min(EthNegative) max(EthPositive)],[z_gnd,z_gnd]*1e-3,'k');
% semilogx(...
%     abs(EzNum(:)), [0:size(EzNum)-1]*dz,'r-',...
%     abs(Einitiation(:)),[0:size(Einitiation)-1]*dz, 'b--',...
%     abs(EthNegative(:)),[0:size(EthNegative)-1]*dz, 'g--',...
%     abs(EthPositive(:)),[0:size(EthPositive)-1]*dz, 'g--'...
%     );

axis([min(EthNegative) max(EthPositive) 0 (Lz+z_gnd)*1e-3]);
legend('E_z [BF]','\phi [BF]','E_{th}^\pm');
xlabel('E_z (V/m) & \phi (kV)');
ylabel('Altitude (km)','FontSize',16);
title('E_z and \phi  after the flash [AF]');
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
box on
hold off


linewidth = 1;
figure;
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(phiNumBF(:)*1e-6, z*1e-3, 'b-','LineWidth',linewidth)
plot(phiNumAF(:)*1e-6, z*1e-3, 'b--','LineWidth',linewidth)
plot([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6)],[z_gnd,z_gnd]*1e-3,'k','LineWidth',linewidth);
hold off
axis([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6) 0 (Lz+z_gnd)*1e-3]);
xlabel('\phi (MV)','FontSize',16);
ylabel('Altitude (km)','FontSize',16);
set(gca,'FontSize',14);
set(gca,'XMinorTick','on','YMinorTick','on')
box on
% title('E_z and \phi  before the flash [BF] and after [AF]');
grid on

figure;
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(EzNumBF(:)*1e-5, z*1e-3, 'r-','LineWidth',linewidth)
plot(EzNumAF(:)*1e-5, z*1e-3, 'r--','LineWidth',linewidth)
plot(EthNegative(:)*1e-5,z*1e-3, 'g-.','LineWidth',linewidth)
plot(EthPositive(:)*1e-5,z*1e-3, 'g-.','LineWidth',linewidth)
plot([1.25*min(EthNegative*1e-5) -1.25*min(EthNegative*1e-5)],[z_gnd,z_gnd]*1e-3,'k','LineWidth',linewidth);
hold off
axis([1.25*min(EthNegative*1e-5) -1.25*min(EthNegative*1e-5) 0 (Lz+z_gnd)*1e-3]);
xlabel('Ez (kV/cm)','FontSize',16);
ylabel('Altitude (km)','FontSize',16);
set(gca,'FontSize',14);
set(gca,'XMinorTick','on','YMinorTick','on')
box on
% title('E_z and \phi  before the flash [BF] and after [AF]');
grid on

fprintf('Average field reduction: %f\n',mean(abs(EzNumAF)./abs(EzNumBF)));
end