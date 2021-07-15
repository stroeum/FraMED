function []=Plot3D_rho()
% Plot an axisymmetric view of the charge layers assumed cylindrical at the
% center of the domain prior to the flash.
clearvars
close all
clc
global Charges NbLayers

%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%

cd ../results
Charges = load('ChargeLayers.dat');
dxyz    = load('dxyz.dat');
Nxyz    = load('Nxyz.dat');
z_gnd   = load('z_gnd.dat');
cd ../viz

gnd_color = [.5 .5 .5];

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = 100*Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3)/100;           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

clear Nxyz
clear dxyz

NbLayers = size(Charges)
NbLayers = NbLayers(1);

%-------------------------------------------------------------------------%
% Draw the layers                                                         %
%-------------------------------------------------------------------------%

z        = (0:Nz-1)'*dz;
disk(Nz) = 0;
for k=1:Nz
    disk(k) = ro(z(k));
end


x             = disk;
y             = z;
[~,rho,z] = cart2pol(x, zeros(size(x,1), size(x,2)), y); 
r             = 25;
theta         = repmat( linspace(0,2*pi,r)', 1, length(rho(:)) );
rho           = repmat(rho(:).', r, 1);
z             = repmat(z(:).', r, 1);
[X,Y,Z]       = pol2cart(theta,rho,z);
% waterfall(X,Y,Z);
hold on
mesh(X+Lx/2,Y+Ly/2,Z+z_gnd,'FaceColor','none');

Xp = [Lx 0 0 Lx]; 
Yp = [Ly Ly 0 0]; 
Zp = [z_gnd z_gnd z_gnd z_gnd];
patch(Xp, Yp, Zp, z_gnd,'FaceColor',gnd_color);
hold off
for n=1:NbLayers
     text(Lx/2+max(Charges(:,5)),Ly/2-max(Charges(:,5)),Charges(n,4)+z_gnd,...
         ['\leftarrow ', num2str(Charges(n,1)),' C'],...
       'HorizontalAlignment','left')
end

axis([0 Lx 0 Ly 0 Lz+z_gnd])
colormap([0 0 0])

end

function [answer]=ro(z)
global Charges NbLayers
% Disk radius %
answer = 0;
for n=1:NbLayers
    if(z>=(Charges(n,4)-Charges(n,6)/2) && z<=(Charges(n,4)+Charges(n,6)/2))
        answer = Charges(n,5);
    end
end
end

% function []=main()
% clear all
% close all
% clc
% 
% ChargeLayers = load('ChargeLayers.dat');
% dxyz         = load('dxyz.dat');
% Nxyz         = load('Nxyz.dat');
% 
% dz = dxyz(3);
% Nz = Nxyz(3);
% z  = [0:Nz-1]'*dz;
% for k=1:Nz
%     disk(k) = ro(z(k));
% end
% 
% 
% x = disk;
% y = z;
% [theta,rho,z] = cart2pol(x, zeros(size(x,1), size(x,2)), y); 
% r = 25;
% theta = repmat( linspace(0,2*pi,r)', 1, length(rho(:)) );
% rho = repmat(rho(:).', r, 1);
% z = repmat(z(:).', r, 1);
% [X,Y,Z] = pol2cart(theta,rho,z);
% % waterfall(X,Y,Z);
% mesh(X+12e3,Y+12e3,Z,'FaceColor','none');
% 
% text(14e+3,10e+3,7e+3,'\leftarrow  31.5 C',...
%      'HorizontalAlignment','left')
% text(14e+3,10e+3,4e+3,'\leftarrow -45.0 C',...
%      'HorizontalAlignment','left')
% text(14e+3,10e+3,1.5e+3,'\leftarrow   4.5 C',...
%      'HorizontalAlignment','left')
% 
% 
% 
% axis([0 24e3 0 24e3 0 12e3])
% colormap([0 0 0])
% 
% end
% 
% function [answer]=ro(z)
% % Disk radius %
% if(z>=(7-1)*1e3 & z<=(7+1)*1e3)
%     answer = 2.5*1e3;
% elseif (z>=(4-1)*1e3 & z<=(4+1)*1e3)
%     answer = 2.0*1e3;
% elseif (z>=(1.5-1)*1e3 & z<=(1.5+1)*1e3)
%     answer = 2.0*1e3;
% else
%     answer = 0*1e3;
% end
% end