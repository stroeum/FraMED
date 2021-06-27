clear all
close all
clc


dxyz         = load('dxyz.dat')*1e-3;
Nxyz         = load('Nxyz.dat');
test         = load('test.dat')*1e-5;
phiNumBF     = load('phiNumBF.dat')*1e-5;
z_gnd        = load('z_gnd.dat')*1e-3;
dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz

plot(phiNumBF,z,'b',test,z,'r--')