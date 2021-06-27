% LMA data loading
% global d Init Links N Nb z_gnd

dxyz              = load('dxyz.dat')*1e-3;                                 %_km
Nxyz              = load('Nxyz.dat');                                      %_dimensionless
InitPoint         = load('InitPoint.dat')*1e-3;                            %_km
Links.data        = load('EstablishedLinks.dat');
z_gnd             = load('z_gnd.dat')*1e-3;                                %_km
Layers.data       = load('ChargeLayers.dat')*1e-3;                         %_kC, _km
rho.YZ            = load('rhoAmbYZ.dat');
rho.XZ            = load('rhoAmbXZ.dat');

% Derive main parameters
Nb.Links          = size(Links.data);
Nb.Links          = Nb.Links(1);
Nb.Points         = Nb.Links+1;
Nb.Layers         = size(Layers.data);
Nb.Layers         = Nb.Layers(1);

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

x                    = (0:N.x-1)*d.x;
y                    = (0:N.y-1)*d.y;
z                    = (0:N.z-1)*d.z+z_gnd;