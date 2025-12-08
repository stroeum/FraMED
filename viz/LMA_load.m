% LMA data loading
% global d Init Links N Nb sims.domain.gnd

cd ../results
InitPoint   = load('InitPoint.dat');                            %_m
Links.data  = load('EstablishedLinks.dat');
Layers.data = load('ChargeLayers.dat');                         %_C, _m
rho.YZ      = load('rhoAmbYZ.dat');
rho.XZ      = load('rhoAmbXZ.dat');
cd ../viz

% Derive main parameters
Nb.Links    = size(Links.data);
Nb.Links    = Nb.Links(1);
Nb.Points   = Nb.Links+1;
Nb.Layers   = size(Layers.data);
Nb.Layers   = Nb.Layers(1);

Init.z      = InitPoint(3);

clear InitPoint

x           = (0:sims.domain.Nx-1)*sims.domain.dx*sims.spatialFactor.Number;
y           = (0:sims.domain.Ny-1)*sims.domain.dy*sims.spatialFactor.Number;
z           = ((0:sims.domain.Nz-1)*sims.domain.dz+sims.domain.gnd)*sims.spatialFactor.Number;