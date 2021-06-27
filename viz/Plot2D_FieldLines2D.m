close all
clear all
clc
% Load datas
close all
clear all
clc

dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
Ex2D         = load('Ex2D.dat');
Ey2D         = load('Ey2D.dat');
Ez2D         = load('Ez2D.dat');
phi2D        = load('Ez2D.dat');
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

fprintf('SOR_data loaded and saved\n')

%% Fix plot parameters

Pcolor         = [1 0 0];
Ncolor         = [0 0 1];
LPcolor        = [1 0 0];
ChargeLayersLineStyle  = '-';
ChargeLayersLineWidth  = 1;

%% Derive data for plotting
x        = (0:Nx-1)'*dx*1e-3;
y        = (0:Ny-1)'*dy*1e-3;
z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
[x,z]    = meshgrid(x,z);

%% Plot
figure;
hold on
% set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/6)
axis([min(x(:)) max(x(:)) min(z(:)) max(z(:))]);

streamslice(x,z,Ex2D',Ez2D',100,'noarrows');
set(findobj('Type','line'),'Color','k')

plot([min(x(:)) max(x(:))], [z_gnd z_gnd]*1e-3,'k');
rectangle('Position',[(ChargeLayers(1,2)-2*ChargeLayers(1,5)/2)*1e-3,(z_gnd+ChargeLayers(1,4)-ChargeLayers(1,7)/2)*1e-3,2*ChargeLayers(1,5)*1e-3,ChargeLayers(1,7)*1e-3],...
    'Curvature',[0,0],...
    'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',LPcolor);
% text(1.05*Lx/2+ChargeLayers(1,5)*1e-3,(ChargeLayers(1,4)+z_gnd)*1e-3,...
%     [num2str(ChargeLayers(1,1),3),' C'],...
%     'HorizontalAlignment','left','BackgroundColor',[1 1 1],...
%     'FontSize',10,'Color',LPcolor)

rectangle('Position',[(ChargeLayers(2,2)-2*ChargeLayers(2,5)/2)*1e-3,(z_gnd+ChargeLayers(2,4)-ChargeLayers(2,7)/2)*1e-3,2*ChargeLayers(2,5)*1e-3,ChargeLayers(2,7)*1e-3],...
    'Curvature',[0,0],...
    'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Ncolor);
% text(1.05*Lx/2+ChargeLayers(2,5)*1e-3,(ChargeLayers(2,4)+z_gnd)*1e-3,...
%     [num2str(ChargeLayers(2,1),3),' C'],...
%     'HorizontalAlignment','left','BackgroundColor',[1 1 1],...
%     'FontSize',10,'Color',Ncolor)

rectangle('Position',[(ChargeLayers(3,2)-2*ChargeLayers(3,5)/2)*1e-3,(z_gnd+ChargeLayers(3,4)-ChargeLayers(3,7)/2)*1e-3,2*ChargeLayers(3,5)*1e-3,ChargeLayers(3,7)*1e-3],...
    'Curvature',[0,0],...
    'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Pcolor);
% text(1.05*Lx/2+ChargeLayers(3,5)*1e-3,(ChargeLayers(3,4)+z_gnd)*1e-3,...
%     [num2str(ChargeLayers(3,1),3),' C'],...
%     'HorizontalAlignment','left','BackgroundColor',[1 1 1],...
%     'FontSize',10,'Color',Pcolor)
hold off
xlabel('x (km)','FontSize',12);
ylabel('z (km)','FontSize',12);
%title('Projection of E-field lines','FontSize',12,'FontWeight','bold');
set(gca,'FontSize',10);
% axis square
box on

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