function [] = Plot2D_FieldLines(num)
cd ../results
dxyz         = load('dxyz.dat');
Nxyz         = load('Nxyz.dat');
Ex2D         = load(['Ex2D',num2str(num),'.dat']);
Ey2D         = load(['Ey2D',num2str(num),'.dat']);
Ez2D         = load(['Ez2D',num2str(num),'.dat']);
phi2D        = load(['phi2D',num2str(num),'.dat']);
z_gnd        = load('z_gnd.dat');
ChargeLayers = load('ChargeLayers.dat');
cd ../viz
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

Color = ['r' 'b' 'r' 'b'];
NbChargeLayers = size(ChargeLayers);
NbChargeLayers = NbChargeLayers(1);
ChargeLayersLineStyle  = '-';
ChargeLayersLineWidth  = 1;

%% Derive data for plotting
x        = (0:Nx-1)'*dx*1e-3;
y        = (0:Ny-1)'*dy*1e-3;
z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
[y,z]    = meshgrid(y,z);

%% Plot
figure;
hold on
%set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/3)
axis([min(y(:)) max(y(:)) min(z(:)) max(z(:))]);


streamslice(y,z,Ey2D,Ez2D,10,'arrows');
% quiver(y,z,Ey2D',Ez2D');
set(findobj('Type','line'),'Color','k')
plot([min(y(:)) max(y(:))], [z_gnd z_gnd]*1e-3,'k');
for ii=1:NbChargeLayers
    rectangle('Position',[(ChargeLayers(ii,3)-2*ChargeLayers(ii,6)/2)*1e-3,(z_gnd+ChargeLayers(ii,4)-ChargeLayers(ii,7)/2)*1e-3,2*ChargeLayers(ii,6)*1e-3,ChargeLayers(ii,7)*1e-3],...
        'Curvature',[0,0],...
        'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',Color(ii));
%     text((ChargeLayers(ii,3)+ChargeLayers(ii,6))*1e-3,(z_gnd+ChargeLayers(ii,4))*1e-3,...
%         ['\leftarrow',num2str(ChargeLayers(ii,1),3),' C'],...
%         'HorizontalAlignment','left','BackgroundColor','w',...
%         'FontSize',10,'Color',Color(ii))
end

hold off
xlabel('y (km)','FontSize',12);
ylabel('z (km)','FontSize',12);
set(gca,'FontSize',10);
axis xy
axis image
box on