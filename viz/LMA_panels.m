% Subplot dimensioning
if (FocusArea.x(1) == 0 && FocusArea.x(2) == 0 && FocusArea.y(1) == 0 && FocusArea.y(2) == 0 && FocusArea.z(1) == 0 && FocusArea.z(2) == 0)
    FocusArea.x(1) = 0;
    FocusArea.x(2) = L.x;
    FocusArea.y(1) = 0;
    FocusArea.y(2) = L.y;
    FocusArea.z(1) = z_gnd;
    FocusArea.z(2) = L.z+z_gnd;
end

% Subplots dimensions %%
A(1:9,:) = [[(FocusArea.y(2) - FocusArea.y(1)) -(FocusArea.x(2) - FocusArea.x(1)) 0 0 0 0 0 0 0 0];...
    [-(FocusArea.z(2) - FocusArea.z(1)) 0 (FocusArea.x(2) - FocusArea.x(1)) 0 0 0 0 0 0 0];...
    [0 -1 0 1 0 0 0 0 0 0];...
    [-1 0 0 0 1 0 0 0 0 0];...
    [0 -(FocusArea.z(2) - FocusArea.z(1)) 0 0 0 (FocusArea.y(2) - FocusArea.y(1)) 0 0 0 0];...
    [0 0 -1 0 0 0 1 0 0 0];...
    [0 0 0 0 0 -1 0 1 0 0];...
    [-1 0 -1 0 0 0 0 0 1 0];...
    [1 0 1 0 0 0 0 0 0 0]];
B(1:9,1) = [0;0;0;0;0;0;0; Figure.hspace(2); Figure.Width - Figure.hspace(1) - Figure.hspace(2) - Figure.hspace(3)];

if (h.tz.position(4) == 0) % by default dimensionate t-z windows to fill the figure frame
    A(10,:) = [0 1 0 1 0 0 0 0 0 1];
    B(10,1) = Figure.Length - Figure.vspace(1) - Figure.vspace(2)- Figure.vspace(3) - Figure.vspace(4);
else                       % use specified dimension
    A(10,:) = [0 0 0 0 0 0 0 0 0 1];
    B(10,1) = h.tz.position(4);
end
S = A\B;

% Assignment of derived dimensions  %
h.xy.position(3)     = S(1);
h.xy.position(4)     = S(2);
h.yz.position(3)     = S(3);
h.yz.position(4)     = S(4);
h.xz.position(3)     = S(5);
h.xz.position(4)     = S(6);
h.nz.position(3)     = S(7);
h.nz.position(4)     = S(8);
h.tz.position(3)     = S(9);
h.tz.position(4)     = S(10);
clear A B S
h.xy.position(1)     = Figure.hspace(1);
h.xy.position(2)     = Figure.vspace(1);
h.yz.position(1)     = h.xy.position(1)+h.xy.position(3)+Figure.hspace(2);
h.yz.position(2)     = h.xy.position(2);
h.xz.position(1)     = h.xy.position(1);
h.xz.position(2)     = h.xy.position(2)+h.xy.position(4)+Figure.vspace(2);
h.nz.position(1)     = h.yz.position(1);
h.nz.position(2)     = h.xz.position(2);
h.tz.position(1)     = h.xy.position(1);
h.tz.position(2)     = h.xz.position(2)+h.xz.position(4)+Figure.vspace(3);

% Plot of empty panels %
h.fig=figure('units','inches','outerposition',[0 0 Figure.Width Figure.Length],'Resize','off');
h.xy.fig = subplot('Position',[0 0 1 1]);
set(h.xy.fig,'Units','inches','Position',h.xy.position);
h.yz.fig = subplot('Position',[1 1 1 .1]);
set(h.yz.fig,'Units','inches','Position',h.yz.position);
h.xz.fig = subplot('Position',[1 1 .1 .1]);
set(h.xz.fig,'Units','inches','Position',h.xz.position);
h.nz.fig = subplot('Position',[1 1 .1 .1]);
set(h.nz.fig,'Units','inches','Position',h.nz.position);
h.tz.fig = subplot('Position',[1 1 .1 .1]);
set(h.tz.fig,'Units','inches','Position',h.tz.position);