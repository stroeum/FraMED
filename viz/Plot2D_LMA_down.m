close all
clear all
clc

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

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%
Nb.Links          = size(Links.data);
Nb.Links          = Nb.Links(1);
Nb.Points         = Nb.Links+1;
Nb.ChargeLayers   = size(ChargeLayers.data);
Nb.ChargeLayers   = Nb.ChargeLayers(1);

Init.z            = InitPoint(3);
N.x               = Nxyz(1);
N.y               = Nxyz(2);
N.z               = Nxyz(3);

d.x               = dxyz(1);                                               %_km
d.y               = dxyz(2);                                               %_km
d.z               = dxyz(3);                                               %_km

L.x               = (N.x-1)*d.x;                                           %_km
L.y               = (N.y-1)*d.y;                                           %_km
L.z               = (N.z-1)*d.z;                                           %_km
clear Nxyz dxyz InitPoint

ChargeLayers.Line.Style = '--';
ChargeLayers.Line.Width = 1;
ChargeLayers.Edge.Color = [[.75 .75 .75];[.75 .75 .75];[.75 .75 .75];[.75 .75 .75]];
% ChargeLayers.Edge.Color = ['none';'none';'none';'none'];
Ground.Line.Width       = .25;

%-------------------------------------------------------------------------%
% Zoom Area                                                               %
%-------------------------------------------------------------------------%
FocusArea.x(1)    = 0;
FocusArea.x(2)    = L.x;
FocusArea.y(1)    = 0;
FocusArea.y(2)    = L.y;
FocusArea.z(1)    = z_gnd;
FocusArea.z(2)    = L.z+z_gnd-2;

%-------------------------------------------------------------------------%
% Fonts                                                                   %
%-------------------------------------------------------------------------%
Links.Line.Width  = 1;
Links.Color       = colormap(jet(Nb.Points));
Links.MarkerSize  = 7;
Font.Size.Labels  = 12;
Font.Size.Axis    = 10;
Font.Name         = 'Helvetica';

%-------------------------------------------------------------------------%
% Altitude histogra                                                       %
%-------------------------------------------------------------------------%
alt_histogram.data(N.z)               = 0;
alt_histogram.data(Links.data(1,3)+1) = 1;
for ii = 1:Nb.Links
    if(Links.data(ii,6)*d.z <= Init.z)
        alt_histogram.data(Links.data(ii,6)+1) = alt_histogram.data(Links.data(ii,6)+1) +1;
    end
end
alt_histogram.max = max(alt_histogram.data);
z                 = ((0:N.z-1)*d.z+z_gnd);
Nb.LowerPoints    = sum(alt_histogram.data);

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
FigureAspectRatio               = 30/38.8235;
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

h.xy.label.x.string   = 'x (km)';
h.xy.label.x.position = [h.xy.position(1)+h.xy.position(3)/2 h.xy.position(2)-(3/4)*Screen.Vertical.Interspace(1) 0 0];
h.yz.label.x.string   = 'z (km)';
h.yz.label.x.position = [h.yz.position(1)+h.yz.position(3)/2 h.xy.label.x.position(2) 0 0];
h.xz.label.x.string   = '';
h.xz.label.x.position = [h.xy.label.x.position(1) h.xz.position(2)-(3/4)*Screen.Vertical.Interspace(2) 0 0];
h.nz.label.x.string   = '';
h.nz.label.x.position = [h.yz.label.x.position(1) h.xz.label.x.position(2) 0 0];
h.tz.label.x.string   = 'step';
h.tz.label.x.position = [h.tz.position(1)+h.tz.position(3)/2 h.tz.position(2)-(3/4)*Screen.Vertical.Interspace(3) 0 0];

h.xy.label.y.string   = 'y (km)';
h.xy.label.y.position = [h.xy.position(1)-(3/5)*Screen.Horizontal.Interspace(1) h.xy.position(2)+h.xy.position(4)/2 0 0];
h.yz.label.y.string   = '';
h.yz.label.y.position = [h.yz.position(1)-(1/2)*Screen.Horizontal.Interspace(2) h.xy.label.y.position(2) 0 0];
h.xz.label.y.string   = 'z (km)';
h.xz.label.y.position = [h.xy.label.y.position(1) h.xz.position(2)+h.xz.position(4)/2 0 0];
h.nz.label.y.string   = '';
h.nz.label.y.position = [h.yz.label.y.position(1) h.xz.label.y.position(2) 0 0];
h.tz.label.y.string   = 'z (km)';
h.tz.label.y.position = [h.xy.label.y.position(1) h.tz.position(2)+h.tz.position(4)/2 0 0];

h.nz.title.string     = '';%'alt-histogram';
h.nz.title.position   = [h.nz.label.x.position(1) h.nz.position(2)+h.nz.position(4)+(1/4)*Screen.Vertical.Interspace(3) 0 0];
h.tz.title.string     = '';%datestr(now, 'yyyymmdd');
h.tz.title.position   = [h.tz.label.x.position(1) h.tz.position(2)+h.tz.position(4)+(1/4)*Screen.Vertical.Interspace(3) 0 0];

h.nz.note.string      = [int2str(Nb.LowerPoints),' pts'];
h.nz.note.position    = [h.nz.position(1)+h.nz.position(3) h.nz.position(2)+h.nz.position(4) 0 0];

%% Subplot bottom left %%
h.xy.fig = subplot('Position',[0 0 1 1]);
set(h.xy.fig,'Units','inches');
set(h.xy.fig,'Position',h.xy.position)

hold on
plot(...
    Links.data(1,1)*d.x,Links.data(1,2)*d.y,...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for ii=1:Nb.Links
    % plots ending point of each link %
    if(Links.data(ii,6)*d.z <= Init.z)
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color',Links.Color(ii,:)...
            );
        axis([FocusArea.x(1) FocusArea.x(2) FocusArea.y(1) FocusArea.y(2)]);
    end
end
for ii=2:2
    rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),2*ChargeLayers.data(ii,5),2*ChargeLayers.data(ii,6)],...
        'Curvature',[1,1],...
        'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
end
hold off

set(h.xy.fig,'FontSize',Font.Size.Axis)
% axis([FocusArea.x(1) FocusArea.x(2) FocusArea.y(1) FocusArea.y(2)]);
box on
annotation(h.fig,'textbox','string',h.xy.label.x.string ,'rotation', 0,...
    'Units','inches','Position',h.xy.label.x.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.xy.label.y.string ,'rotation',90,...
    'Units','inches','Position',h.xy.label.y.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');

%% Subplot bottom right %%
h.yz.fig = subplot('Position',[1 1 1 .1]);
set(h.yz.fig,'Units','inches');
set(h.yz.fig,'Position',h.yz.position)

hold on
plot(...
    (Links.data(1,3)*d.z+z_gnd),Links.data(1,2)*d.y,...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for ii=1:Nb.Links
    % plots ending point of each link %
    if(Links.data(ii,6)*d.z <= Init.z)
        plot(...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            [Links.data(ii,2)*d.y, Links.data(ii,5)*d.y],...
            'Color',Links.Color(ii,:)...
            );
        axis([FocusArea.z(1) FocusArea.z(2) FocusArea.y(1) FocusArea.y(2)]);
    end
end
plot([z_gnd, z_gnd], [FocusArea.y(1) FocusArea.y(2)], 'k', 'LineWidth', Ground.Line.Width);
for ii=1:Nb.ChargeLayers-1
    rectangle('Position',[(z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)/2),(ChargeLayers.data(ii,3)-ChargeLayers.data(ii,6)),ChargeLayers.data(ii,7),2*ChargeLayers.data(ii,6)],...
        'Curvature',[0,0],...
        'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
end
hold off

% set(h.yz.fig,'YTickLabel','');
set(h.yz.fig,'FontSize',Font.Size.Axis);
% axis([FocusArea.z(1) FocusArea.z(2) FocusArea.y(1) FocusArea.y(2)]);
box on
annotation(h.fig,'textbox','string',h.yz.label.x.string ,'rotation', 0,...
    'Units','inches','Position',h.yz.label.x.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.yz.label.y.string ,'rotation',90,...
    'Units','inches','Position',h.yz.label.y.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');

%% Subplot middle left %%
h.xz.fig = subplot('Position',[1 1 .1 .1]);
set(h.xz.fig,'Units','inches');
set(h.xz.fig,'Position',h.xz.position)

hold on
plot(...
    Links.data(1,1)*d.x,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for ii=1:Nb.Links
    % plots ending point of each link %
    if(Links.data(ii,6)*d.z <= Init.z)
        plot(...
            [Links.data(ii,1)*d.x, Links.data(ii,4)*d.x],...
            [Links.data(ii,3)*d.z+z_gnd, Links.data(ii,6)*d.z+z_gnd],...
            'Color',Links.Color(ii,:)...
            );
        axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);
    end
end
plot([FocusArea.x(1) FocusArea.x(2)], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);
for ii=1:Nb.ChargeLayers-1
    rectangle('Position',[(ChargeLayers.data(ii,2)-ChargeLayers.data(ii,5)),...
        (z_gnd+ChargeLayers.data(ii,4)-ChargeLayers.data(ii,7)/2),....
        2*ChargeLayers.data(ii,5),ChargeLayers.data(ii,7)],...
        'Curvature',[0,0],...
        'LineWidth',ChargeLayers.Line.Width,'LineStyle',ChargeLayers.Line.Style,'EdgeColor',ChargeLayers.Edge.Color(ii,:));
end
hold off

% set(h.xz.fig,'XTickLabel','');
set(h.xz.fig,'FontSize',Font.Size.Axis);
% axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);
box on
annotation(h.fig,'textbox','string',h.xz.label.x.string ,'rotation', 0,...
    'Units','inches','Position',h.xz.label.x.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.xz.label.y.string ,'rotation',90,...
    'Units','inches','Position',h.xz.label.y.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');

%% Subplot middle right %%
h.nz.fig = subplot('Position',[1 1 .1 .1]);
set(h.nz.fig,'Units','inches');
set(h.nz.fig,'Position',h.nz.position)

plot(alt_histogram.data,z,'k-');
axis([0 alt_histogram.max FocusArea.z(1) FocusArea.z(2)]);

% set(h.nz.fig,'XTickLabel','');
set(h.nz.fig,'FontSize',Font.Size.Axis);
box on
annotation(h.fig,'textbox','string',h.nz.label.x.string ,'rotation', 0,...
    'Units','inches','Position',h.nz.label.x.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.nz.label.y.string ,'rotation',90,...
    'Units','inches','Position',h.nz.label.y.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.nz.title.string ,'rotation', 0,...
    'Units','inches','Position',h.nz.title.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.nz.note.string ,'rotation', 0,...
    'Units','inches','Position',h.nz.note.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','top','HorizontalAlignment','right');

%% Subplot top %%
h.tz.fig = subplot('Position',[1 1 .1 .1]);
set(h.tz.fig,'Units','inches');
set(h.tz.fig,'Position',h.tz.position)

hold on
plot(...
    0,(Links.data(1,3)*d.z+z_gnd),...
    'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,...
    'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));
for n=1:Nb.Links
    % plots ending point of each link %
    if(Links.data(n,6)*d.z <= Init.z)
        plot(...
            n,(Links.data(n,6)*d.z+z_gnd),...
            's-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),...
            'MarkerEdgeColor',Links.Color(n,:),'MarkerFaceColor',Links.Color(n,:));
        axis([0 Nb.Points FocusArea.z(1) FocusArea.z(2)]);
    end
end
hold off

set(h.tz.fig,'FontSize',Font.Size.Axis)
% axis([0 Nb.Points FocusArea.z(1) FocusArea.z(2)]);
box on
annotation(h.fig,'textbox','string',h.tz.label.x.string ,'rotation', 0,...
    'Units','inches','Position',h.tz.label.x.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.tz.label.y.string ,'rotation',90,...
    'Units','inches','Position',h.tz.label.y.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');
annotation(h.fig,'textbox','string',h.tz.title.string ,'rotation', 0,...
    'Units','inches','Position',h.tz.title.position,...
    'EdgeColor','none','FontSize',Font.Size.Labels,'FontName',Font.Name,...
    'VerticalAlignment','bottom','HorizontalAlignment','center');

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