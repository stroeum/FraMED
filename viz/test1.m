clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
EzNumBF      = load('EnumBF.dat');
phiNumBF     = load('phiNumBF.dat');
EzNumAF      = load('EnumAF.dat');
phiNumAF     = load('phiNumAF.dat');
% Einitiation  = load('Einitiation.dat');
EthPositive  = load('EthPositive.dat');
EthNegative  = load('EthNegative.dat');
z_gnd        = load('z_gnd.dat');
phiMonopole  = load('Monopole.dat');
phiDipole    = load('Dipole.dat');
phiMultipole = load('Multipole.dat');
dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linewidth = 1;
figure;
subplot(121)
% set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(log10(phiNumBF(:)*1e-6),     z*1e-3, 'r-', 'LineWidth',linewidth)
plot(log10(phiMonopole(:)*1e-6),  z*1e-3, 'g--','LineWidth',linewidth)
plot(log10(phiDipole(:)*1e-6),    z*1e-3, 'b--','LineWidth',linewidth)
plot(log10(phiMultipole(:)*1e-6), z*1e-3, 'k--','LineWidth',linewidth)
legend('Numerical','Monopole','Dipole','Multipole');
% plot([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6)],[z_gnd,z_gnd]*1e-3,'k','LineWidth',linewidth);
hold off
% axis([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6) 0 (Lz+z_gnd)*1e-3]);
xlabel('\phi (MV)','FontSize',16);
ylabel('Altitude (km)','FontSize',16);
set(gca,'FontSize',14);
box on
% title('E_z and \phi  before the flash [BF] and after [AF]');
grid on

subplot(122)
% set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
eta1 = 100*abs(phiNumBF-phiMonopole)  ./ phiMonopole;
eta2 = 100*abs(phiNumBF-phiDipole)    ./ phiDipole;
eta3 = 100*abs(phiNumBF-phiMultipole) ./ phiMultipole;
plot(eta1, z*1e-3, 'g-', 'LineWidth',linewidth)
plot(eta2, z*1e-3, 'b-', 'LineWidth',linewidth)
plot(eta3, z*1e-3, 'k-', 'LineWidth',linewidth)
axis([0 max(eta1(2:end-1)) min(z*1e-3) max(z*1e-3)]);
legend('Monopole','Dipole','Multipole');
hold off
xlabel('$\eta [\%] = \frac{|\phi_{num}-\phi_{ana}|}{\phi_{ana}}$','Interpreter','Latex','FontSize',16);
ylabel('Altitude (km)','FontSize',16);
set(gca,'FontSize',14);
box on
% title('E_z and \phi  before the flash [BF] and after [AF]');
grid on
