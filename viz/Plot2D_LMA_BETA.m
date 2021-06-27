close all
clear all
clc
beep  off

%% Load relevant data %%
%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%

dxyz              = load('dxyz.dat')*1e-3;                                 %_km
Nxyz              = load('Nxyz.dat');                                      %_dimensionless
InitPoint         = load('InitPoint.dat')*1e-3;                          %_km
Links.data        = load('EstablishedLinks.dat');
z_gnd             = load('z_gnd.dat')*1e-3;                                %_km
ChargeLayers.data = load('ChargeLayers.dat')*1e-3;                         %_kC, _km
rhoYZ             = load('rhoAmbYZ.dat');
rhoXZ             = load('rhoAmbXZ.dat');

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%
Nb.Links          = size(Links.data);
Nb.Links          = Nb.Links(1);
Nb.Points         = Nb.Links+1;
Nb.ChargeLayers   = size(ChargeLayers.data);
Nb.ChargeLayers   = Nb.ChargeLayers(1);

N.x               = Nxyz(1);
N.y               = Nxyz(2);
N.z               = Nxyz(3);

d.x               = dxyz(1);                                               %_km
d.y               = dxyz(2);                                               %_km
d.z               = dxyz(3);                                               %_km

Init.z            = InitPoint(3);
L.x               = (N.x-1)*d.x;                                           %_km
L.y               = (N.y-1)*d.y;                                           %_km
L.z               = (N.z-1)*d.z;                                           %_km
clear Nxyz dxyz InitPoint

ChargeLayers.Type       = 'none'; % 'disk'; % 'spheres'; %
ChargeLayers.Line.Style = '--';
ChargeLayers.Line.Width = 1;
ChargeLayers.Edge.Color = [[.75 .75 .75];[.75 .75 .75];[.75 .75 .75];[.75 .75 .75]];
% ChargeLayers.Edge.Color = ['none';'none';'none';'none'];
ChargeLayers.Density    = 'on';
Ground.Line.Width       = .25;

%-------------------------------------------------------------------------%
% Zoom Area                                                               %
%-------------------------------------------------------------------------%
FocusArea.x(1)    = 0;
FocusArea.x(2)    = L.x;
FocusArea.y(1)    = 0;
FocusArea.y(2)    = L.y;
FocusArea.z(1)    = z_gnd;
FocusArea.z(2)    = L.z+z_gnd;

%-------------------------------------------------------------------------%
% Fonts                                                                   %
%-------------------------------------------------------------------------%
Links.Line.Width  = 1;
Links.Scheme      = 'Colormap';
Links.MarkerSize  = 7;
Font.Size.Labels  = 12;
Font.Size.Axis    = 10;
Font.Name         = 'Helvetica';

if strcmp(Links.Scheme,'Colormap')
    Links.Color       = colormap(jet(Nb.Points));
elseif strcmp(Links.Scheme,'B/W')
    Links.Color       = [0 0 0];
elseif strcmp(Links.Scheme,'R/B')
    Links.Color(1,:)  = [0 0 0];
    Links.Color(2,:)  = [1 0 0];
    Links.Color(3,:)  = [0 0 1];
elseif strcmp(Links.Scheme,'B/R')
    Links.Color(1,:)  = [0 0 0];
    Links.Color(2,:)  = [0 0 1];
    Links.Color(3,:)  = [1 0 0];
end

%-------------------------------------------------------------------------%
% Altitude histogra                                                       %
%-------------------------------------------------------------------------%
alt_histogram.data(N.z)               = 0;
alt_histogram.up(N.z)                 = 0;
alt_histogram.down(N.z)               = 0;
alt_histogram.data(Links.data(1,3)+1) = 1;
if(Links.data(1,3)*d.z >= Init.z)
    alt_histogram.up(Links.data(1,3)+1)   = 1;
elseif (Links.data(1,3)*d.z <= Init.z)
    alt_histogram.down(Links.data(1,3)+1) = 1;
end
for ii = 1:Nb.Links
    alt_histogram.data(Links.data(ii,6)+1) = alt_histogram.data(Links.data(ii,6)+1) +1;
    if(Links.data(ii,6)*d.z >= Init.z)
        alt_histogram.up(Links.data(ii,6)+1)   = alt_histogram.up(Links.data(ii,6)+1) +1;
    elseif (Links.data(ii,6)*d.z <= Init.z)
        alt_histogram.down(Links.data(ii,6)+1) = alt_histogram.down(Links.data(ii,6)+1) +1;
    end
end
alt_histogram.max = max(alt_histogram.data);
x                 = (0:N.x-1)*d.x;
y                 = (0:N.y-1)*d.y;
z                 = (0:N.z-1)*d.z+z_gnd;

%% Preliminary dimensioning %%
h.fig = figure('units','normalized','outerposition',[0 0 1 1]);            % Create full screen figure
set(h.fig,'Units','inch');                                                 % Set units used for the formatting of the figure to Inches
H1     = get(h.fig,'Outerposition');
H2     = get(h.fig,'Position');
HL     = H2(4);
dH     = H1(4)-H2(4);
close all
clear h H1 H2

%% Resize figure dimension %%
FigureAspectRatio               = .6;%30/38.8235;
Paper.Figure.Width              = 30;                                      %_picas
Paper.Horizontal.Interspace(1)  = .511811;                                 %_inches
Paper.Horizontal.Interspace(2)  = .275591;                                 %_inches
Paper.Horizontal.Interspace(3)  = .157480;                                 %_inches
Paper.Vertical.Interspace(1)    = .472441;                                 %_inches
Paper.Vertical.Interspace(2)    = .236220;                                 %_inches
Paper.Vertical.Interspace(3)    = .236220;                                 %_inches
Paper.Vertical.Interspace(4)    = .433071-.275591-.078740;                 %_inches

Paper.Figure.Length             = Paper.Figure.Width/FigureAspectRatio;    %_picas
Paper.Figure.Width              = Paper.Figure.Width/6;                    %_inches
Paper.Figure.Length             = Paper.Figure.Length/6;                   %_inches

h.fig=figure('units','normalized','outerposition',[0 0 1 .9]);             % Create full screen figure
set(h.fig,'Units','inch');                                                 % Set units used for the formatting of the figure to Inches
PrintableArea                   = get(h.fig,'Outerposition');              % Get maximum printable size in inches
Screen.Figure.Length            = PrintableArea(4);                        %_inches: Length of the figure on the screen
Screen.Figure.Width             = Screen.Figure.Length*FigureAspectRatio;  %_inches: Length of the figure on the screen

set(h.fig,'Outerposition',[0 0 Screen.Figure.Width Screen.Figure.Length]); % redimensionate figure frame to maximize figure dimension on the screen

ScreenPaper.LengthRatio         = Screen.Figure.Length/Paper.Figure.Length;
ScreenPaper.WidthRatio          = Screen.Figure.Width/Paper.Figure.Width;
ScreenPaper.Ratio               = ScreenPaper.LengthRatio;                 % Ratio between the size of the figure on the screen and the size of the figure on paper

%% Subplots dimensions %%
if(ScreenPaper.Ratio <=1)
    Screen.Horizontal.Interspace(1) = Paper.Horizontal.Interspace(1)*ScreenPaper.Ratio; %_inches
    Screen.Horizontal.Interspace(2) = Paper.Horizontal.Interspace(2)*ScreenPaper.Ratio; %_inches
    Screen.Horizontal.Interspace(3) = Paper.Horizontal.Interspace(3)*ScreenPaper.Ratio; %_inches
    Screen.Vertical.Interspace(1)   = Paper.Vertical.Interspace(1)*ScreenPaper.Ratio;   %_inches
    Screen.Vertical.Interspace(2)   = Paper.Vertical.Interspace(2)*ScreenPaper.Ratio;   %_inches
    Screen.Vertical.Interspace(3)   = Paper.Vertical.Interspace(3)*ScreenPaper.Ratio;   %_inches
    Screen.Vertical.Interspace(4)   = Paper.Vertical.Interspace(4)*ScreenPaper.Ratio;   %_inches
else
    Screen.Horizontal.Interspace(1) = Paper.Horizontal.Interspace(1);                   %_inches
    Screen.Horizontal.Interspace(2) = Paper.Horizontal.Interspace(2);                   %_inches
    Screen.Horizontal.Interspace(3) = Paper.Horizontal.Interspace(3);                   %_inches
    Screen.Vertical.Interspace(1)   = Paper.Vertical.Interspace(1);                     %_inches
    Screen.Vertical.Interspace(2)   = Paper.Vertical.Interspace(2);                     %_inches
    Screen.Vertical.Interspace(3)   = Paper.Vertical.Interspace(3);                     %_inches
    Screen.Vertical.Interspace(4)   = Paper.Vertical.Interspace(4);                     %_inches
    Screen.Figure.Width             = Paper.Figure.Width;                               %_inches: Length of the figure on the screen
    Screen.Figure.Length            = Screen.Figure.Width/FigureAspectRatio;            %_inches: Length of the figure on the screen
    set(h.fig,'Outerposition',[0 0 Screen.Figure.Width Screen.Figure.Length]);
end
h.xy.position                   = zeros(1,4);
h.yz.position                   = zeros(1,4);
h.xz.position                   = zeros(1,4);
h.nz.position                   = zeros(1,4);
h.tz.position                   = zeros(1,4);

h.xy.label.x.string             = 'x (km)';
h.xy.label.y.string             = 'y (km)';
h.xy.label.spacing.x            = 16;
h.xy.label.spacing.y            = 24;
h.yz.label.x.string             = 'z (km)';
h.yz.label.y.string             = '';%'y (km)';
h.yz.label.spacing.x            = h.xy.label.spacing.x;
h.yz.label.spacing.y            = h.xy.label.spacing.y;
h.yz.label.spacing.x            = h.xy.label.spacing.x;
h.yz.label.spacing.y            = h.xy.label.spacing.y;
h.xz.label.x.string             = '';%'x (km)';
h.xz.label.y.string             = 'z (km)';
h.xz.label.spacing.x            = h.xy.label.spacing.x;
h.xz.label.spacing.y            = h.xy.label.spacing.y;
h.nz.label.spacing.x            = h.xy.label.spacing.x;
h.nz.label.x.string             = '';%'n';
h.nz.label.y.string             = '';%'z (km)';
h.nz.title.string               = '';%'alt-histogram';
h.nz.note.string                = [int2str(Nb.Points),' pts '];
h.nz.label.spacing.y            = h.xy.label.spacing.y;
h.tz.label.x.string             = 'step';
h.tz.label.y.string             = 'z (km)';
h.tz.title.string               = '';%datestr(now, 'yyyymmdd');
h.tz.label.spacing.x            = h.xy.label.spacing.x;
h.tz.label.spacing.y            = h.xy.label.spacing.y;

B = [0;0;0;0;0;0;0;Screen.Horizontal.Interspace(2);...
    Screen.Figure.Width - Screen.Horizontal.Interspace(1) - Screen.Horizontal.Interspace(2) ...
    - Screen.Horizontal.Interspace(3);...
    Screen.Figure.Length - Screen.Vertical.Interspace(1) - Screen.Vertical.Interspace(2) ...
    - Screen.Vertical.Interspace(3) - Screen.Vertical.Interspace(4)];
A = [[(FocusArea.y(2) - FocusArea.y(1)) -(FocusArea.x(2) - FocusArea.x(1)) 0 0 0 0 0 0 0 0];...
    [-(FocusArea.z(2) - FocusArea.z(1)) 0 (FocusArea.x(2) - FocusArea.x(1)) 0 0 0 0 0 0 0];...
    [0 -1 0 1 0 0 0 0 0 0];...
    [-1 0 0 0 1 0 0 0 0 0];...
    [0 -(FocusArea.z(2) - FocusArea.z(1)) 0 0 0 (FocusArea.y(2) - FocusArea.y(1)) 0 0 0 0];...
    [0 0 -1 0 0 0 1 0 0 0];...
    [0 0 0 0 0 -1 0 1 0 0];...
    [-1 0 -1 0 0 0 0 0 1 0];...
    [1 0 1 0 0 0 0 0 0 0];...
    [0 1 0 1 0 0 0 0 0 1]];

S = A\B;
h.xy.position(3)      = S(1);
h.xy.position(4)      = S(2);
h.yz.position(3)      = S(3);
h.yz.position(4)      = S(4);
h.xz.position(3)      = S(5);
h.xz.position(4)      = S(6);
h.nz.position(3)      = S(7);
h.nz.position(4)      = S(8);
h.tz.position(3)      = S(9);
h.tz.position(4)      = S(10);
clear A B S
h.xy.position(1)      = Screen.Horizontal.Interspace(1);
h.xy.position(2)      = Screen.Vertical.Interspace(1);
h.yz.position(1)      = h.xy.position(1)+h.xy.position(3)+Screen.Horizontal.Interspace(2);
h.yz.position(2)      = h.xy.position(2);
h.xz.position(1)      = h.xy.position(1);
h.xz.position(2)      = h.xy.position(2)+h.xy.position(4)+Screen.Vertical.Interspace(2);
h.nz.position(1)      = h.yz.position(1);
h.nz.position(2)      = h.xz.position(2);
h.tz.position(1)      = h.xy.position(1);
h.tz.position(2)      = h.xz.position(2)+h.xz.position(4)+Screen.Vertical.Interspace(3);

%% Subplot bottom left %%
h.xy.fig = subplot('Position',[0 0 1 1]);
set(h.xy.fig,'Units','inches');
set(h.xy.fig,'Position',h.xy.position)
set(h.xy.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(h.xy.fig,'XMinorTick','on','YMinorTick','on')

hold on
plot(...
    Links.data(1,1)*d.x,Links.data(1,2)*d.y,...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for ii=1:Nb.Links
    % plots ending point of each link %
    if strcmp(Links.Scheme,'Colormap')
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color',Links.Color(ii,:)...
            );
    elseif strcmp(Links.Scheme,'B/W')
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color','k'...
            );
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if(Links.data(ii,6)*d.z >= Init.z)
            plot(...
                [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                'Color',Links.Color(2,:)...
                );
        else
            plot(...
                [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                'Color',Links.Color(3,:)...
                );
        end
    end
    axis([FocusArea.x(1) FocusArea.x(2) FocusArea.y(1) FocusArea.y(2)]);
end
if      strcmp(ChargeLayers.Type,'disks')
    for ii=1:Nb.ChargeLayers
        rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),2*ChargeLayers.data(ii,5),2*ChargeLayers.data(ii,6)],...
            'Curvature',[1,1],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
elseif  strcmp(ChargeLayers.Type,'spheres')
    for ii=2:3
        rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),2*ChargeLayers.data(ii,5),2*ChargeLayers.data(ii,6)],...
            'Curvature',[1,1],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
end
hold off
% axis([FocusArea.x(1) FocusArea.x(2) FocusArea.y(1) FocusArea.y(2)]);
box on

h.xy.label.x.textbox     = xlabel(h.xy.label.x.string,'Units','Points');
h.xy.label.x.position    = get(h.xy.label.x.textbox,'Position');
h.xy.label.x.position(2) = -h.xy.label.spacing.x;
set(h.xy.label.x.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.xy.label.x.position)
h.xy.label.y.textbox     = ylabel(h.xy.label.y.string,'Units','Points');
h.xy.label.y.position    = get(h.xy.label.y.textbox,'Position');
h.xy.label.y.position(1) = -h.xy.label.spacing.y;
set(h.xy.label.y.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.xy.label.y.position)

%% Subplot bottom right %%
h.yz.fig = subplot('Position',[1 1 1 .1]);
set(h.yz.fig,'Units','inches');
set(h.yz.fig,'Position',h.yz.position)
set(h.yz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(h.yz.fig,'XMinorTick','on','YMinorTick','on')
% set(h.yz.fig,'YTickLabel','');

hold on
if strcmp(ChargeLayers.Density,'on')
    contourf(z,y,rhoYZ,60,'LineColor','none');
    caxis([-max(max(abs(rhoYZ))) max(max(abs(rhoYZ)))])
    % imagesc(Y,Z,rhoYZ');
%     colorbar('Location','East');
    box on
    set(gca,'XMinorTick','on','YMinorTick','on')
end
plot(...
    (Links.data(1,3)*d.z+z_gnd),Links.data(1,2)*d.y,...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for ii=1:Nb.Links
    % plots ending point of each link %
    if strcmp(Links.Scheme,'Colormap')
        plot(...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color',Links.Color(ii,:)...
            );
    elseif strcmp(Links.Scheme,'B/W')
        plot(...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color','k'...
            );
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if(Links.data(ii,6)*d.z >= Init.z)
            plot(...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                'Color',Links.Color(2,:)...
                );
        else
            plot(...
                [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
                [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
                'Color',Links.Color(3,:)...
                );
        end
    end
    axis([FocusArea.z(1) FocusArea.z(2) FocusArea.y(1) FocusArea.y(2)]);
end
plot([z_gnd, z_gnd], [FocusArea.y(1) FocusArea.y(2)], 'k', 'LineWidth', Ground.Line.Width);
if      strcmp(ChargeLayers.Type,'disks')
    for ii=1:Nb.ChargeLayers
        rectangle('Position',[(z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)/2),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),ChargeLayers.data(ii,7),2*ChargeLayers.data(ii,6)],...
            'Curvature',[0,0],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
elseif  strcmp(ChargeLayers.Type,'spheres')
    for ii=2:3
        rectangle('Position',[(z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),2*ChargeLayers.data(ii,7),2*ChargeLayers.data(ii,6)],...
            'Curvature',[1,1],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
end
hold off

% axis([FocusArea.z(1) FocusArea.z(2) FocusArea.y(1) FocusArea.y(2)]);
box on
h.yz.label.x.textbox     = xlabel(h.yz.label.x.string,'Units','Points');
h.yz.label.x.position    = get(h.yz.label.x.textbox,'Position');
h.yz.label.x.position(2) = -h.yz.label.spacing.x;
set(h.yz.label.x.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.yz.label.x.position)
h.yz.label.y.textbox     = ylabel(h.yz.label.y.string,'Units','Points');
h.yz.label.y.position    = get(h.yz.label.y.textbox,'Position');
h.yz.label.y.position(1) = -h.yz.label.spacing.y;
set(h.yz.label.y.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.yz.label.y.position)

%% Subplot middle left %%
h.xz.fig = subplot('Position',[1 1 .1 .1]);
set(h.xz.fig,'Units','inches');
set(h.xz.fig,'Position',h.xz.position)
set(h.xz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(h.xz.fig,'XMinorTick','on','YMinorTick','on')
% set(h.xz.fig,'XTickLabel','');

hold on
if strcmp(ChargeLayers.Density,'on')
    contourf(x,z,rhoXZ',60,'LineColor','none');
    caxis([-max(max(abs(rhoXZ))) max(max(abs(rhoXZ)))])
    % imagesc(Y,Z,rhoYZ');
%     colorbar('Location','North');
end
plot(...
    Links.data(1,1)*d.x,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
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
    axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);
end
plot([FocusArea.x(1) FocusArea.x(2)], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);

if      strcmp(ChargeLayers.Type,'disks')
    for ii=1:Nb.ChargeLayers
        rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),...
            (z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)/2),....
            2*ChargeLayers.data(ii,5),ChargeLayers.data(ii,7)],...
            'Curvature',[0,0],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
elseif  strcmp(ChargeLayers.Type,'spheres')
    for ii=2:3
        rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),...
            (z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)),....
            2*ChargeLayers.data(ii,5),2*ChargeLayers.data(ii,7)],...
            'Curvature',[1,1],...
            'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
    end
end
hold off

% axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);
box on
h.xz.label.x.textbox     = xlabel(h.xz.label.x.string,'Units','Points');
h.xz.label.x.position    = get(h.xz.label.x.textbox,'Position');
h.xz.label.x.position(2) = -h.xz.label.spacing.x;
set(h.xz.label.x.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.xz.label.x.position)
h.xz.label.y.textbox     = ylabel(h.xz.label.y.string,'Units','Points');
h.xz.label.y.position    = get(h.xz.label.y.textbox,'Position');
h.xz.label.y.position(1) = -h.xz.label.spacing.y;
set(h.xz.label.y.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.xz.label.y.position)

%% Subplot middle right %%
h.nz.fig = subplot('Position',[1 1 .1 .1]);
set(h.nz.fig,'Units','inches');
set(h.nz.fig,'Position',h.nz.position)
set(h.nz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(h.nz.fig,'XMinorTick','on','YMinorTick','on')

if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
    plot(alt_histogram.data,z,'k-');
elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
    hold on
    plot(alt_histogram.up,z,'Color',Links.Color(2,:));
    plot(alt_histogram.down,z,'Color',Links.Color(3,:));
    plot(alt_histogram.data,z,'k--');
    hold off
end
axis([0 alt_histogram.max FocusArea.z(1) FocusArea.z(2)]);

box on
h.nz.label.x.textbox     = xlabel(h.nz.label.x.string,'Units','Points');
h.nz.label.x.position    = get(h.nz.label.x.textbox,'Position');
h.nz.label.x.position(2) = -h.nz.label.spacing.x;
set(h.nz.label.x.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.nz.label.x.position)
h.nz.label.y.textbox     = ylabel(h.nz.label.y.string,'Units','Points');
h.nz.label.y.position    = get(h.nz.label.y.textbox,'Position');
h.nz.label.y.position(1) = -h.nz.label.spacing.y;
set(h.nz.label.y.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.nz.label.y.position)
h.nz.title.textbox       = title(h.nz.title.string,'Units','normalized');
h.nz.title.position      = get(h.nz.title.textbox,'Position');
h.nz.title.position(2)   = 1;
set(h.nz.title.textbox,...
    'Units','normalized','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Bottom','Position',h.nz.title.position)
text(alt_histogram.max,FocusArea.z(2),h.nz.note.string,...
    'FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Top','HorizontalAlignment','Right')

%% Subplot top %%
h.tz.fig = subplot('Position',[1 1 .1 .1]);
set(h.tz.fig,'Units','inches');
set(h.tz.fig,'Position',h.tz.position)
set(h.tz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name)
set(h.tz.fig,'XMinorTick','on','YMinorTick','on')

hold on
plot(...
    0,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for n=1:Nb.Links
    % plots ending point of each link %
    if strcmp(Links.Scheme,'Colormap')
        plot(...
            n,(Links.data(n,6)*d.z+z_gnd),...
            's-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),...
            'MarkerEdgeColor',Links.Color(n,:),'MarkerFaceColor',Links.Color(n,:));
    elseif strcmp(Links.Scheme,'B/W')
        plot(...
            n,(Links.data(n,6)*d.z+z_gnd),...
            's-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if(Links.data(n,6)*d.z >= Init.z)
            plot(...
                n,(Links.data(n,6)*d.z+z_gnd),...
                's-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),...
                'MarkerEdgeColor',Links.Color(2,:),'MarkerFaceColor',Links.Color(2,:));
        else
            plot(...
                n,(Links.data(n,6)*d.z+z_gnd),...
                's-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),...
                'MarkerEdgeColor',Links.Color(3,:),'MarkerFaceColor',Links.Color(3,:));
        end
    end
    axis([0 Nb.Points FocusArea.z(1) FocusArea.z(2)]);
end
hold off

% axis([0 Nb.Points FocusArea.z(1) FocusArea.z(2)]);
box on
h.tz.label.x.textbox     = xlabel(h.tz.label.x.string,'Units','Points');
h.tz.label.x.position    = get(h.tz.label.x.textbox,'Position');
h.tz.label.x.position(2) = 0;
set(h.tz.label.x.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Top','Position',h.tz.label.x.position)
h.tz.label.y.textbox     = ylabel(h.tz.label.y.string,'Units','Points');
h.tz.label.y.position    = get(h.tz.label.y.textbox,'Position');
h.tz.label.y.position(1) = -h.tz.label.spacing.y;
set(h.tz.label.y.textbox,...
    'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Middle','Position',h.tz.label.y.position)
h.tz.title.textbox       = title(h.tz.title.string,'Units','normalized');
h.tz.title.position      = get(h.tz.title.textbox,'Position');
h.tz.title.position(2)   = 1;
set(h.tz.title.textbox,...
    'Units','normalized','FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Bottom','Position',h.tz.title.position)

%% Final Resizing %%
H = get(h.fig,'Outerposition');
HW = H(3);
if(ScreenPaper.Ratio >=1)
    HL = H(4)+.75*dH;                                                       % if PaperSize of the figure is bigger than the screen, then use the max size of the screen, else redimensionate adequately
end
set(h.fig,'Outerposition',[0,0,HW, HL])
clear H HW HL ii n z

% get(h.fig,'OuterPosition')
% get(h.xy.fig,'Position')
% get(h.yz.fig,'Position')
% get(h.xz.fig,'Position')
% get(h.tz.fig,'Position')
% get(h.nz.fig,'Position')
% set(gcf,'HandleVisibility','off')