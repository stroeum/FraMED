function []=Plot1D_rho()
clear all
close all
clc
global Charges NbLayers

cd ../results/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
Charges     = load('ChargeLayers.dat');
z_gnd       = load('z_gnd.dat');

dz = dxyz(3)/100;
Nz = Nxyz(3)*100;
Lz = (Nz-1)*dz;
NbLayers = size(Charges);
NbLayers = NbLayers(1)-1;

clear Nxyz dxyz

z       = (0:Nz-1)'*dz;
rho(Nz) = 0;
for k=1:Nz
    rho(k) = ro(z(k));
end

hold on
plot(rho,z+z_gnd,'r');
plot([0 max(rho)],[z_gnd,z_gnd],'k');
hold off

xlabel('Charge extension');
ylabel('Altitude (m)');
grid on
axis([0 max(rho) 0 Lz+z_gnd])

cd ../viz/
end

function [answer]=ro(z)
global Charges NbLayers
% Disk radius %
answer = 0;
for n=1:NbLayers
    if(z>=(Charges(n,4)-Charges(n,6)/2) && z<=(Charges(n,4)+Charges(n,6)/2))
        answer = Charges(n,5);
    end
end
end