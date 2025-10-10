% LMA data loading
% global d Init Links N Nb z_gnd

cd ../results
dxyz              = load('dxyz.dat');                                 %_m
Nxyz              = load('Nxyz.dat');                                 %_dimensionless
InitPoint         = load('InitPoint.dat');                            %_m
Links.data        = load('EstablishedLinks.dat');
z_gnd             = load('z_gnd.dat');                                %_m
Layers.data       = load('ChargeLayers.dat');                         %_C, _m
rho.YZ            = load('rhoAmbYZ.dat');
rho.XZ            = load('rhoAmbXZ.dat');
cd ../viz

% Determine the magnitude of the values:
spatialFactor = checkMagnitude((0:(Nxyz(3)-1))*dxyz(3)+z_gnd);

% Derive main parameters
Nb.Links          = size(Links.data);
Nb.Links          = Nb.Links(1);
Nb.Points         = Nb.Links+1;
Nb.Layers         = size(Layers.data);
Nb.Layers         = Nb.Layers(1);

N.x               = Nxyz(1);
N.y               = Nxyz(2);
N.z               = Nxyz(3);

d.x               = dxyz(1)*spatialFactor.Number;                     %_converted unit
d.y               = dxyz(2)*spatialFactor.Number;                     %_converted unit
d.z               = dxyz(3)*spatialFactor.Number;                     %_converted unit
z_gnd             = z_gnd*spatialFactor.Number;                       %_converted unit

Init.z            = InitPoint(3);
L.x               = (N.x-1)*d.x;                                      %_converted unit
L.y               = (N.y-1)*d.y;                                      %_converted unit
L.z               = (N.z-1)*d.z;                                      %_converted unit

clear Nxyz dxyz InitPoint

x                    = (0:N.x-1)*d.x;
y                    = (0:N.y-1)*d.y;
z                    = (0:N.z-1)*d.z+z_gnd;