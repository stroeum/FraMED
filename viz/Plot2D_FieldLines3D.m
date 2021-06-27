close all
clear all
clc
%% Load datas
close all
clear all
clc

dxyz         = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/dxyz.dat');
Nxyz         = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/Nxyz.dat');
Ex           = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/Ex.dat');
Ey           = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/Ey.dat');
Ez           = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/Ez.dat');
z_gnd        = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/z_gnd.dat');
ChargeLayers = load('/Users/Jeremy/Documents/Academic/Research/2. Code/FractalModel_Random_Init/results/ChargeLayers.dat');

dx = dxyz(1);               % _m
dy = dxyz(2);               % _m
dz = dxyz(3);               % _m

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

Ex  = ConvertTo3d(Ex ,Nxyz);
Ey  = ConvertTo3d(Ey ,Nxyz);
Ez  = ConvertTo3d(Ez ,Nxyz);

clear dxyz Nxyz

save('SOR_Efield.mat')
fprintf('SOR_data loaded and saved\n')

% load SOR_Efield

%% Fix plot parameters

% LineWidth      = 1;
% MarkerSize     = 10;
% FontSizeLabels = 16;
% FontSizeAxis   = 12;
Pcolor         = [1 0 0];
Ncolor         = [0 0 1];
LPcolor        = [1 0 0];
ChargeLayersLineStyle  = '-';
ChargeLayersLineWidth  = 1;

%% Derive data for plotting
x = (0:Nx-1)'*dx*1e-3;
y = (0:Ny-1)'*dy*1e-3;
z = ((0:Nz-1)'*dz+z_gnd)*1e-3;

u(Nx,Nz) = 0; % x-component of E-field

v(Nx,Nz) = 0; % h-component of E-field (horizontal)
w(Nx,Nz) = 0; % v-component of E-field (vertical)
for k = 1:Nz
    for i = 1:Nx
        u(i,k) = Ex(i,(Ny+1)/2,k);
        if(i>(Nx+1)/2)
            v(i,k) = sqrt(Ex(i,(Ny+1)/2,k)^2+Ey(i,(Ny+1)/2,k)^2);
        else
            v(i,k) = -sqrt(Ex(i,(Ny+1)/2,k)^2+Ey(i,(Ny+1)/2,k)^2);
        end
        w(i,k) = Ez(i,(Ny+1)/2,k);
    end
end
[x,z]=meshgrid(x,z);

%% Plot
figure;
axis([min(x(:)) max(x(:)) 0 max(z(:))]);
hold on
streamslice(x,z,u',w','Color','k');
plot([min(x(:)) max(x(:))], [z_gnd z_gnd]*1e-3,'k');
rectangle('Position',[(ChargeLayers(1,2)-2*ChargeLayers(1,5)/2)*1e-3,(z_gnd+ChargeLayers(1,4)-ChargeLayers(1,6)/2)*1e-3,2*ChargeLayers(1,5)*1e-3,ChargeLayers(1,6)*1e-3],...
          'Curvature',[0,0],...
         'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',LPcolor);
rectangle('Position',[(ChargeLayers(2,2)-2*ChargeLayers(2,5)/2)*1e-3,(z_gnd+ChargeLayers(2,4)-ChargeLayers(2,6)/2)*1e-3,2*ChargeLayers(2,5)*1e-3,ChargeLayers(2,6)*1e-3],...
          'Curvature',[0,0],...
         'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Ncolor);
rectangle('Position',[(ChargeLayers(3,2)-2*ChargeLayers(3,5)/2)*1e-3,(z_gnd+ChargeLayers(3,4)-ChargeLayers(3,6)/2)*1e-3,2*ChargeLayers(3,5)*1e-3,ChargeLayers(3,6)*1e-3],...
         'Curvature',[0,0],...
         'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Pcolor);
hold off

% figure;
% hold on
% streamslice(x,z,v',w','Color','k');
% rectangle('Position',[(ChargeLayers(1,2)-2*ChargeLayers(1,5)/2)*1e-3,(z_gnd+ChargeLayers(1,4)-ChargeLayers(1,6)/2)*1e-3,2*ChargeLayers(1,5)*1e-3,ChargeLayers(1,6)*1e-3],...
%           'Curvature',[0,0],...
%          'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',LPcolor);
% rectangle('Position',[(ChargeLayers(2,2)-2*ChargeLayers(2,5)/2)*1e-3,(z_gnd+ChargeLayers(2,4)-ChargeLayers(2,6)/2)*1e-3,2*ChargeLayers(2,5)*1e-3,ChargeLayers(2,6)*1e-3],...
%           'Curvature',[0,0],...
%          'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Ncolor);
% rectangle('Position',[(ChargeLayers(3,2)-2*ChargeLayers(3,5)/2)*1e-3,(z_gnd+ChargeLayers(3,4)-ChargeLayers(3,6)/2)*1e-3,2*ChargeLayers(3,5)*1e-3,ChargeLayers(3,6)*1e-3],...
%          'Curvature',[0,0],...
%          'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Pcolor);
% hold off