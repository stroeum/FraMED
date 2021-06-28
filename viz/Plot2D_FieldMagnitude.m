%% Load datas
close all
clear all
clc

%% Compile
% cd ..
% !make
% !./main
% cd results

%% Draw
dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
Ex2D         = load('Ex2D.dat');
Ey2D         = load('Ey2D.dat');
Ez2D         = load('Ez2D.dat');
Einitiation  = load('Einitiation.dat');
phi2D        = load('phi2D.dat');
z_gnd        = load('z_gnd.dat');
ChargeLayers = load('ChargeLayers.dat');

dx = dxyz(1);               % _m
dy = dxyz(2);               % _m
dz = dxyz(3);               % _m

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

Lx = (Nx-1)*dx*1e-3;         % _m
Ly = (Ny-1)*dy*1e-3;         % _m
Lz = (Nz-1)*dz*1e-3;         % _m

clear dxyz Nxyz

fprintf('Data loaded\n')

%% Fix plot parameters

NbChargeLayers = size(ChargeLayers);
NbChargeLayers = NbChargeLayers(1);
ChargeLayersLineStyle  = '-';
ChargeLayersLineWidth  = 1;

%% Derive data for plotting
x        = (0:Nx-1)'*dx*1e-3;
y        = (0:Ny-1)'*dy*1e-3;
z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
% [x,z]    = meshgrid(x,z);

E        = (Ex2D.^2+Ey2D.^2+Ez2D.^2).^.5*1e-5;
[M, N]   = size(E);
Einit    = zeros(M,N);
for mm=1:M
    Einit(mm,:) = Einitiation(:)*1e-5;
end

%% Plot

% set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/3)

hold on
isFieldMagnitude = 2; %input('Plot E-field magnitude? (1: yes, else: no)\n>> ');
if (isFieldMagnitude ~= 1)
    imagesc(y,z,(((E-Einit)./Einit)'));
    title('Threshold Overshoot','FontSize',12,'FontWeight','bold');
elseif (isFieldMagnitude == 1)
    imagesc(y,z,E');
    title('E-field magnitude (kV/cm)','FontSize',12,'FontWeight','bold');
end
colormap(gray(256))
colorbar

% for ii=1:NbChargeLayers
%     rectangle('Position',[(ChargeLayers(ii,3)-2*ChargeLayers(ii,6)/2)*1e-3,(z_gnd+ChargeLayers(ii,4)-ChargeLayers(ii,7)/2)*1e-3,2*ChargeLayers(ii,6)*1e-3,ChargeLayers(ii,7)*1e-3],...
%         'Curvature',[0,0],...
%         'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor','w');
% %     text((ChargeLayers(ii,3)+ChargeLayers(ii,6))*1e-3,(z_gnd+ChargeLayers(ii,4))*1e-3,...
% %         ['\leftarrow',num2str(ChargeLayers(ii,1),3),' C'],...
% %         'HorizontalAlignment','left','BackgroundColor','none',...
% %         'FontSize',10,'Color','w')
% end
xlabel('y (km)','FontSize',12);
ylabel('z (km)','FontSize',12);
set(gca,'FontSize',10);
axis xy
axis image
% axis([0 Ly 0 25]);
box on
hold off