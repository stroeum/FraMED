close all
clearvars -except sims


cd ../results/
figure;
%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%

load dxyz.dat               -ascii
load Nxyz.dat               -ascii
load InitPoint.dat          -ascii
load EstablishedLinks.dat   -ascii
load z_gnd.dat              -ascii
load ChargeLayers.dat       -ascii

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%
NbOfLinks = size(EstablishedLinks);
NbOfLinks = NbOfLinks(1);

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m
LzMin = z_gnd;
LzMax = (Lz+z_gnd);
% InitX = InitPoint(1);   % _m
% InitY = InitPoint(2);   % _m
% InitZ = InitPoint(3);   % _m
% InitR = InitPoint(4);   % _m

% Initi = round(InitX/dx);
% Initj = round(InitY/dy);
% Initk = round(InitZ/dz);

clear Nxyz
clear dxyz
clear InitPoint

X = (0:Nx-1)*dx;
Y = (0:Ny-1)*dy;
% Z = (0:Nz-1)*dz;

NbOfPoints              = NbOfLinks+1;
% color                   = colormap(jet(NbOfPoints));
% color                   = colormap(gray(1.25*NbOfPoints));
color     = colormap(jet(NbOfLinks));
isMonochrome = 1;%input('Is plot Monochrome? (1: yes, else: no)\n>> ');
if (isMonochrome == 1)
    for ii=1:NbOfLinks
        color(ii,:) = [0 0 1];
    end
end
LineWidth               = 1;
MarkerSize              = 10;
FontSizeLabels          = 16;
FontSizeAxis            = 12;
Pcolor                  = [.75 .75 .75];
Ncolor                  = [.75 .75 .75];
LPcolor                 = [.75 .75 .75];
ChargeLayersLineStyle  = '-';
ChargeLayersLineWidth  = 1;

%-------------------------------------------------------------------------%
% Derive subplots parameters                                              %
%-------------------------------------------------------------------------%

%% Derive data for plotting
x        = (0:Nx-1)'*dx*1e-3;
y        = (0:Ny-1)'*dy*1e-3;
z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
[y,z]    = meshgrid(y,z);

%% Plot
hold on
%set(gcf,'Units','inches','OuterPosition', [10 10 20 20]/6)
set(gcf,'Units','inches','Position', [10 10 20 10]/6)

% axis([min(y(:)) max(y(:)) min(z(:)) max(z(:))]);

%set(gca,'Units','pixels','FontSize',FontSizeAxis)
box on
hold on;

rectangle('Position',[(ChargeLayers(1,3)-2*ChargeLayers(1,5)/2)*1e-3,(z_gnd+ChargeLayers(1,4)-ChargeLayers(1,6)/2)*1e-3,2*ChargeLayers(1,5)*1e-3,ChargeLayers(1,6)*1e-3],...
    'Curvature',[0,0],...
    'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',LPcolor);
rectangle('Position',[(ChargeLayers(2,3)-2*ChargeLayers(2,5)/2)*1e-3,(z_gnd+ChargeLayers(2,4)-ChargeLayers(2,6)/2)*1e-3,2*ChargeLayers(2,5)*1e-3,ChargeLayers(2,6)*1e-3],...
    'Curvature',[0,0],...
    'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Ncolor);
% rectangle('Position',[(ChargeLayers(3,3)-2*ChargeLayers(3,5)/2)*1e-3,(z_gnd+ChargeLayers(3,4)-ChargeLayers(3,6)/2)*1e-3,2*ChargeLayers(3,5)*1e-3,ChargeLayers(3,6)*1e-3],...
%           'Curvature',[0,0],...
%          'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Pcolor);

plot(Y*1e-3,z_gnd*1e-3*ones(size(Y)),'k','LineWidth',LineWidth);
plot(...
    EstablishedLinks(1,2)*dy*1e-3,(EstablishedLinks(1,3)*dz+z_gnd)*1e-3,...
    'x-','LineWidth',LineWidth,'MarkerSize',MarkerSize,...
    'MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:));
for ii=1:NbOfLinks
    % plots ending point of each link %
    plot(...
        [EstablishedLinks(ii,2)*dy, EstablishedLinks(ii,5)*dy]*1e-3,...
        [EstablishedLinks(ii,3)*dz+z_gnd, EstablishedLinks(ii,6)*dz+z_gnd]*1e-3,...
        'Color',color(ii,:)...
        );
    axis([0 Ly LzMin LzMax]*1e-3);
    axis equal tight xy
    %      axis([0 40 0 20]);
end
xlabel('y (km)','FontSize',FontSizeLabels)
ylabel('z (km)','FontSize',FontSizeLabels)

hold off;

cd ../viz/
