close all
clearvars -except sims

beep  off

%% Load relevant data %%
%-------------------------------------------------------------------------%
% Load data files %
%-------------------------------------------------------------------------%
cd ../results/
dxyz = load('dxyz.dat')*1e-3; %_km
Nxyz = load('Nxyz.dat'); %_dimensionless
InitPoint = load('InitPoint.dat')*1e-3; %_km
Links.data = load('EstablishedLinks.dat');
z_gnd = load('z_gnd.dat')*1e-3; %_km
ChargeLayers.data = load('ChargeLayers.dat')*1e-3; %_kC, _km
rho.YZ = load('rhoAmbYZ.dat');
rho.XZ = load('rhoAmbXZ.dat');

%-------------------------------------------------------------------------%
% Derive main parameters %
%-------------------------------------------------------------------------%
Nb.Links = size(Links.data);
Nb.Links = Nb.Links(1);
Nb.Points = Nb.Links+1;
Nb.ChargeLayers = size(ChargeLayers.data);
Nb.ChargeLayers = Nb.ChargeLayers(1);

N.x = Nxyz(1);
N.y = Nxyz(2);
N.z = Nxyz(3);

d.x = dxyz(1); %_km
d.y = dxyz(2); %_km
d.z = dxyz(3); %_km

Init.z = InitPoint(3);
L.x = (N.x-1)*d.x; %_km
L.y = (N.y-1)*d.y; %_km
L.z = (N.z-1)*d.z; %_km
clear Nxyz dxyz InitPoint

ChargeLayers.Type = 'disks' ; %'none' ; % 'spheres'; %
ChargeLayers.Line.Style = '-';
ChargeLayers.Line.Width = 1;
ChargeLayers.FontSize = 8;
ChargeLayers.Pos.x = 10;
ChargeLayers.Circles = 12;
ChargeLayers.Edge.Color = [[1 0 0];[0 0 1];[1 0 0];[.75 .75 .75]];
% ChargeLayers.Edge.Color = ['none';'none';'none';'none'];
ChargeLayers.Density = 'on';
Ground.Line.Width = .25;

%-------------------------------------------------------------------------%
% Zoom Area %
%-------------------------------------------------------------------------%
FocusArea.x(1) = 0;
FocusArea.x(2) = L.x;
FocusArea.y(1) = 0;
FocusArea.y(2) = L.y;
FocusArea.z(1) = z_gnd;
FocusArea.z(2) = L.z+z_gnd;

%-------------------------------------------------------------------------%
% Fonts %
%-------------------------------------------------------------------------%
Links.Line.Width = 1;
Links.Scheme = 'R/B';
Links.MarkerSize = 3;
Font.Size.Labels = 12;
Font.Size.Axis = 10;
Font.Name = 'Helvetica';

if strcmp(Links.Scheme,'Colormap')
    Links.Color = colormap(jet(Nb.Points));
elseif strcmp(Links.Scheme,'B/W')
    Links.Color = [0 0 0];
elseif strcmp(Links.Scheme,'R/B')
    Links.Color(1,:) = [0 0 0];
    Links.Color(2,:) = [1 0 0];
    Links.Color(3,:) = [0 0 1];
elseif strcmp(Links.Scheme,'B/R')
    Links.Color(1,:) = [0 0 0];
    Links.Color(2,:) = [0 0 1];
    Links.Color(3,:) = [1 0 0];
end

x = (0:N.x-1)*d.x;
y = (0:N.y-1)*d.y;
z = (0:N.z-1)*d.z+z_gnd;

%% Plot
figure(1);
set(gcf,'Units','normalized','OuterPosition', [0 1/3 1 2/3])

%% Subplot xz %%
subplot(121);
set(gca,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XAxisLocation','bottom','YAxisLocation','left')
set(gca,'TickDir','in')
set(gca,'TickLength',[.025 .05])
set(gca,'LineWidth',.25)
% set(gca,'YTickLabel','')

hold on
if strcmp(ChargeLayers.Density,'on')
    imagescSgnLog(x,z,rho.XZ',-3,1);
    axis([0 12 z(1)*1e-3 z(end)*1e-3])
    xlabel('x (km)','fontsize',18)
    ylabel('z (km)','fontsize',18)
    XTick = get(gca,'XTick');
    set(gca,'XMinorTick','on','YMinorTick','on');
    title('\rho_t=\rho_s+\rho_f (nC/m^3)','fontsize',18);
    set(gca,'FontSize',18);
end
plot(...
    Links.data(1,1)*d.x,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width/5,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

for ii=1:Nb.Links
    % plots ending point of each link %
    if strcmp(Links.Scheme,'Colormap')
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            'Color',Links.Color(ii,:)...
            );
    elseif strcmp(Links.Scheme,'B/W')
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            'Color','k'...
            );
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if(Links.data(ii,6)*d.z >= Init.z)
            plot(...
                [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                'Color',Links.Color(2,:)...
                );
        else
            plot(...
                [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                'Color',Links.Color(3,:)...
                );
        end
    end
end
plot([FocusArea.x(1) FocusArea.x(2)], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);
hold off

box on
axis image
axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);

%% Subplot yz %%
subplot(122);
set(gca,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'XAxisLocation','bottom','YAxisLocation','left')
set(gca,'TickDir','in')
set(gca,'TickLength',[.025 .05])
set(gca,'LineWidth',.25)
% set(gca,'YTickLabel','')

hold on
if strcmp(ChargeLayers.Density,'on')
    imagescSgnLog(y,z,rho.YZ',-3,1);
    axis([y(1) y(end) z(1)*1e-3 z(end)*1e-3])
    xlabel('y (km)','fontsize',18)
    ylabel('z (km)','fontsize',18)
    XTick = get(gca,'XTick');
    set(gca,'XMinorTick','on','YMinorTick','on');
    title('\rho_t=\rho_s+\rho_f (nC/m^3)','fontsize',18);
    set(gca,'FontSize',18);
end
plot(...
    Links.data(1,2)*d.y,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width/5,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

for ii=1:Nb.Links
    % plots ending point of each link %
    if strcmp(Links.Scheme,'Colormap')
        plot(...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            'Color',Links.Color(ii,:)...
            );
    elseif strcmp(Links.Scheme,'B/W')
        plot(...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            'Color','k'...
            );
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if(Links.data(ii,6)*d.z >= Init.z)
            plot(...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                'Color',Links.Color(2,:)...
                );
        else
            plot(...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                'Color',Links.Color(3,:)...
                );
        end
    end
end
plot([FocusArea.y(1) FocusArea.y(2)], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);
hold off

box on
axis image
axis([FocusArea.y(1) FocusArea.y(2) FocusArea.z(1) FocusArea.z(2)]);
cd ../viz

function [M]=imagescSgnLog(varargin)
% output = imagescSgnLog(density matrix)
% output = imagescSgnLog(r,z,density matrix)
% output = imagescSgnLog(r,z,density matrix,lower boundary)
% output = imagescSgnLog(r,z,density matrix,lower boundary, upper boundary)
% output = imagescSgnLog(r,z,density matrix,lower boundary, upper boundary,colorbar label)
% Plot a logscaled image of the input matrix accounting for the original
% sign of each element of the input
format short e
narginchk(1,6);
if(nargin==2)
    error('Wrong number of arguments (1, 3-6)');
end

if nargin==1
    NN  = size(varargin{1});
    N.r = NN(1);
    N.z = NN(2);
    M.data = varargin{1};
    clear NN;
else
    N.r = length(varargin{1});
    N.z = length(varargin{2});
    M.data = varargin{3};
end

M.log.P = log10((M.data>=0).*M.data);
M.log.N = log10(-(M.data<=0).*M.data);

M.max.P = max(max(M.log.P));
M.max.N = max(max(M.log.N));
M.max.T = max(M.max.P,M.max.N);

M.min.P = min(min(M.log.P));
M.min.N = min(min(M.log.N));
M.min.T = min(M.min.P,M.min.N);

if nargin<4
    if(isinf(M.max.T))
        M.max.T = 2;
    end
    if(isinf(M.min.T))
        M.min.T = -6;
    end
elseif nargin == 4
    M.min.T = varargin{4};
elseif nargin > 4
    M.min.T = varargin{4};
    M.max.T = varargin{5};
end
for ii=1:N.r
    for kk=1:N.z
        if M.log.P(ii,kk) < M.min.T
            M.log.P(ii,kk) = M.min.T;
        end
        if M.log.N(ii,kk) < M.min.T
            M.log.N(ii,kk) = M.min.T;
        end
        
    end
end

M.len.T = ceil(2*(M.max.T-M.min.T));

% M.log.P = M.log.P+M.len.T/2;
% M.log.N = -M.log.N-M.len.T/2;
M.log.P = M.log.P-M.min.T;
M.log.N = -M.log.N+M.min.T;

M.log.T = (M.data>=0).*M.log.P + (M.data<0).*M.log.N;

% subplot(131); imagesc(M.log.P'); colorbar; title('+'); axis xy; axis image; caxis([-M.len.T/2 M.len.T/2]);
% subplot(132); imagesc(M.log.N'); colorbar; title('-'); axis xy; axis image; caxis([-M.len.T/2 M.len.T/2]);
% subplot(133); imagesc(M.log.T'); colorbar; title('T'); axis xy; axis image; caxis([-M.len.T/2 M.len.T/2]);

if nargin==1
    imagesc(M.log.T');
else
    imagesc(varargin{1},varargin{2},M.log.T');
end
axis xy;
axis image;
caxis([-M.len.T/2 M.len.T/2]);
cbar = colorbar;
% cbar.AxisLocation = 'in';
cTick = cbar.Ticks;
c = '';
for ii=1:numel(cTick)
    if cTick(ii)<0
        c{ii} = num2str(-10^(-cTick(ii)+M.min.T),'%2.1e\n');
    elseif cTick(ii)==0
        c{ii} = ['+/-',num2str(10^M.min.T,'%2.1e\n')];
    elseif cTick(ii)>0
        c{ii} = num2str( 10^(cTick(ii)+M.min.T),'%2.1e\n');
    end
    
end
if nargin<6
    cbar.TickLabels = c;
    cbar.TickDirection = 'out';
else
    cbar.TickLabels = c;
    cbar.TickDirection = 'out';
    cbar.Label.String = varargin{6};
end    
end
