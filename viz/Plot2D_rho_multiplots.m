close all
clear all
clc
choice = 0;

dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
z_gnd       = load('z_gnd.dat');
dy = dxyz(2);
dz = dxyz(3);
Ny = Nxyz(2);
Nz = Nxyz(3);
clear Nxyz dxyz


Y = (0:dy:(Ny-1)*dy);
Z = (0:dz:(Nz-1)*dz)+z_gnd;

if choice == 0
    rho2D       = load('rho2d.dat');
    imagesc(Z,Y,rho2D,[-10 10]);
    ylabel('y-axis (m)');
    xlabel('z-axis (m)');
    title(['\rho  (nC) after ', int2str(390) ,' step(s)']);
    colormap(jet(256))
    colorbar;
elseif choice == 1
    M = 3;
    N = 2;
    P = M*N;
    figure;
    SavingStep = 40;
    for k=1:P
        subplot(M,N,k)
        fname=['rho2d',int2str(k*SavingStep),'.dat'];
        rho2D       = load(fname);
        imagesc(Z,Y,rho2D,[-10 10]);
        ylabel('y-axis (m)');
        xlabel('z-axis (m)');
        title(['\rho  (nC) after ', int2str(k*SavingStep) ,' step(s)']);
        colormap(jet(256))
        colorbar;
    end

    figure;
    for k=P+1:12
        subplot(M,N,k-P)
        fname=['rho2d',int2str(k*SavingStep),'.dat'];
        rho2D       = load(fname);
        imagesc(Z,Y,rho2D,[-10 10]);
        ylabel('y-axis (m)');
        xlabel('z-axis (m)');
        title(['\rho  (nC) after ', int2str(k*SavingStep) ,' step(s)']);
        colormap(jet(256))
        colorbar;
    end
end