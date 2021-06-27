% LMA main %
close all
clear all
clc
beep  off
global d Init Links N Nb z_gnd
global alt_hist FocusArea Font Ground

%% Initialisation %%
% !rm -rf *.avi *.eps
Nb.Links = 0; Nb.Points = 0; Nb.Layers = 0;
     d.x = 0;       d.y = 0;       d.z = 0;
  Init.x = 0;    Init.y = 0;    Init.z = 0;
     N.x = 0;       N.y = 0;       N.z = 0;
   z_gnd = 0;

%% Preferences %%
% Links %
Figure.Output     =  'Plot'; % 'Movie'; %

Links.Type        = 'Both'; %Lower'; % 'Upper'; %
Links.Scheme      = 'B/R';
Links.Line.Width  = 1;
Links.MarkerSize  = 7;
Font.Size.Labels  = 12;
Font.Size.Axis    = 10;
Font.Name         = 'Helvetica';

% Charge centers %
Layers.Type       = 'arbitrary'; % 'disks'; % 'spheres'; % 'none';
Layers.Line.Style = '--';
Layers.Line.Width = 1;
Layers.Edge.Color = [[1 0 0];[0 0 1];[1 0 0];[0 0 1]];

% Ground plane %
Ground.Line.Width = .25;

% Zoom area %
FocusArea.x(1)    = 0;
FocusArea.x(2)    = 0;
FocusArea.y(1)    = 0;
FocusArea.y(2)    = 0;
FocusArea.z(1)    = 0;
FocusArea.z(2)    = 0;

% Figure dimension %
Figure.Width       = 5.0;                                                  %_inches
Figure.Length      = 7.5;                                                  %_inches
Figure.hspace(1)   = .511811;                                              %_inches
Figure.hspace(2)   = .275591;                                              %_inches
Figure.hspace(3)   = .157480;                                              %_inches
Figure.vspace(1)   = .472441;                                              %_inches
Figure.vspace(2)   = .236220;                                              %_inches
Figure.vspace(3)   = .236220;                                              %_inches
Figure.vspace(4)   = 1.5;                                                  %_inches
h.tz.position(4)   = 1.0;                                                  %_inches

% Label properties %
h.xy.label.x.string  = 'x (km)';
h.xy.label.y.string  = 'y (km)';
h.xy.label.spacing.x = 16;
h.xy.label.spacing.y = 24;
h.yz.label.x.string  = 'z (km)';
h.yz.label.y.string  = '';%'y (km)';
h.yz.label.spacing.x = h.xy.label.spacing.x;
h.yz.label.spacing.y = h.xy.label.spacing.y;
h.yz.label.spacing.x = h.xy.label.spacing.x;
h.yz.label.spacing.y = h.xy.label.spacing.y;
h.xz.label.x.string  = '';%'x (km)';
h.xz.label.y.string  = 'z (km)';
h.xz.label.spacing.x = h.xy.label.spacing.x;
h.xz.label.spacing.y = h.xy.label.spacing.y;
h.nz.label.spacing.x = h.xy.label.spacing.x;
h.nz.label.x.string  = '';%'n';
h.nz.label.y.string  = '';%'z (km)';
h.nz.title.string    = '';%'alt-histogram';
h.nz.note.string     = [int2str(0),' pts '];
h.nz.label.spacing.y = h.xy.label.spacing.y;
h.tz.label.x.string  = 'step';
h.tz.label.y.string  = 'z (km)';
h.tz.title.string    = '';%datestr(now, 'yyyymmdd');
h.tz.label.spacing.x = h.xy.label.spacing.x;
h.tz.label.spacing.y = h.xy.label.spacing.y;

%% Load parameters %%
LMA_load

%% Create color scheme %%
LMA_color

%% Create figure panels %%
LMA_panels

%% Create altitude histogram %%
alt_hist = LMA_alt_hist(alt_hist,-1);
subplot(h.nz.fig)
h.nz.note.string = [int2str(0),' pts'];
text(alt_hist.max,FocusArea.z(2),h.nz.note.string,'FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Top','HorizontalAlignment','Right')

%% Modify subplot axes %%
LMA_axes(h,'all');

%% Add subplots labels %%
LMA_labels(h,'all')

%% Plot charge regions %%
LMA_charges

%% Plot ground plane %%
LMA_ground(h,'all')

%% Plot initiation point %%
LMA_initiation

%% Plot discharge tree branches %%
LMA_tree

%% Output results %%
if strcmp(Figure.Output,'Movie')
    movie2avi(Movie,'LMA.avi','quality',100);
elseif strcmp(Figure.Output,'Plot')
    hgexport(h.fig,'LMA.eps');
end