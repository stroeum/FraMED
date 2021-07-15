close all
clearvars
clc
choice = 1;
cd ../results/
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
    rho2D       = load('rhoAmbYZ0.dat');
    imagesc(Y*1e-3,Z*1e-3,rho2D,[-10 10]);
    xlabel('y-axis (m)');
    ylabel('z-axis (m)');
    title(['\rho  (nC) after ', int2str(390) ,' step(s)']);
    axis xy equal
    colormap(jet(256))
    colorbar;
elseif choice == 1
    M = 3;
    N = 2;
    P = M*N;
    figure;
    SavingStep = 10;
    for k=1:P
        subplot(M,N,k)
        fname=['rhoAmbYZ',int2str(k*SavingStep),'.dat'];
        rho2D       = load(fname);
        imagesc(Y*1e-3,Z*1e-3,rho2D,[-10 10]);
        ylabel('y-axis (m)');
        xlabel('z-axis (m)');
        title(['\rho  (nC) after ', int2str(k*SavingStep) ,' step(s)']);
        axis xy equal tight
        col 

    figure;
    for k=P+1:12
        subplot(M,N,k-P)
        fname=['rhoAmbYZ',int2str(k*SavingStep),'.dat'];
        rho2D       = load(fname);
        imagesc(Y*1e-3,Z*1e-3,rho2D,[-10 10]);
        ylabel('y-axis (m)');
        xlabel('z-axis (m)');
        title(['\rho  (nC) after ', int2str(k*SavingStep) ,' step(s)']);
        axis xy equal tight
        colormap(jet(256))
        colorbar;
    end
end

cd ../viz/