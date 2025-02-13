function []=Plot3D_TreeAndField()
close all
clearvars -except sims


figure;
%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%

cd ../results
% load Un.dat                 -ascii
load dxyz.dat               -ascii
load Nxyz.dat               -ascii
load InitPoint.dat          -ascii
load EstablishedLinks.dat   -ascii
% EzNum       = load('EzNum.dat');
phiNum      = load('phiNumBF.dat');
% Einitiation = load('Einitiation.dat');
% EthPositive = load('EthPositive.dat');
% EthNegative = load('EthNegative.dat');
z_gnd       = load('z_gnd.dat');

gnd_color = [.5 .5 .5];

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%
NbOfLinks = size(EstablishedLinks);
NbOfLinks = NbOfLinks(1);

Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

InitX = InitPoint(1);   % _m
InitY = InitPoint(2);   % _m
InitZ = InitPoint(3);   % _m
InitR = InitPoint(4);   % _m

Initi = round(InitX/dx);
Initj = round(InitY/dy);
Initk = round(InitZ/dz);

clear Nxyz
clear dxyz
clear InitPoint

%-------------------------------------------------------------------------%
% Draw the tree                                                           %
%-------------------------------------------------------------------------%
figure(1); 
hold on;

%-----------------------------%
% Superimpose initial E-field %
%-----------------------------%
z = (0:Nz-1)*dz;
% x = EzNum/max(max(max(abs(EzNum))))*5e3;
x = phiNum/max(max(max(abs(phiNum))))*5e3;
y = z;
[theta,rho,z] = cart2pol(x, zeros(size(x,1), size(x,2)), y); 
r = 17;
theta = repmat( linspace(0,2*pi,r)', 1, length(rho(:)) );
rho = repmat(rho(:).', r, 1);
z = repmat(z(:).', r, 1);
[X,Y,Z] = pol2cart(theta,rho,z);
% waterfall(X,Y,Z);
mesh(X+Lx/2,Y+Ly/2,Z+z_gnd,'FaceColor','none');
% legend('5 E_{axis} [kV/m]')
legend('5 \phi_{axis} [kV]')
colormap([.75 .75 .75])

% xc = InitX;
% yc = InitY;
% zc = InitZ;
% lx = 4*dx;%0.0025;%2e-3;
% ly = 4*dy;%0.0025;%2e-3;
% lz = 4*dz;%0.0025;%2e-3;
% 
% x1 = xc-lx/2;
% x2 = xc+lx/2;
% y1 = yc-ly/2;
% y2 = yc+ly/2;
% z1 = zc-lz/2;
% z2 = zc+lz/2;
% 
% plot3([x1 x2],[y1 y1],[z1 z1],'k-');
% plot3([x2 x2],[y1 y2],[z1 z1],'k-');
% plot3([x2 x1],[y2 y2],[z1 z1],'k-');
% plot3([x1 x1],[y2 y1],[z1 z1],'k-');
% plot3([x1 x2],[y1 y1],[z2 z2],'k-');
% plot3([x2 x2],[y1 y2],[z2 z2],'k-');
% plot3([x2 x1],[y2 y2],[z2 z2],'k-');
% plot3([x1 x1],[y2 y1],[z2 z2],'k-');
% plot3([x1 x1],[y1 y1],[z1 z2],'k-');
% plot3([x2 x2],[y1 y1],[z1 z2],'k-');
% plot3([x2 x2],[y2 y2],[z1 z2],'k-');
% plot3([x1 x1],[y2 y2],[z1 z2],'k-');
% 
% plot3([x1 x2],[y1 y2],[z1 z1],'k-');
% plot3([x2 x1],[y1 y2],[z1 z1],'k-');
% plot3([x1 x2],[y1 y2],[z2 z2],'k-');
% plot3([x2 x1],[y1 y2],[z2 z2],'k-');
% 
% plot3([x1 x2],[y1 y1],[z1 z2],'k-');
% plot3([x2 x1],[y1 y1],[z1 z2],'k-');
% plot3([x1 x2],[y2 y2],[z1 z2],'k-');
% plot3([x2 x1],[y2 y2],[z1 z2],'k-');
% 
% plot3([x1 x1],[y1 y2],[z1 z2],'k-');
% plot3([x1 x1],[y2 y1],[z1 z2],'k-');
% plot3([x2 x2],[y1 y2],[z1 z2],'k-');
% plot3([x2 x2],[y2 y1],[z1 z2],'k-');
% 
% plot3([x1 x2],[y1 y2],[z1 z2],'k-');
% plot3([x1 x2],[y1 y2],[z2 z1],'k-');
% plot3([x2 x1],[y1 y2],[z1 z2],'k-');
% plot3([x1 x1],[y2 y1],[z1 z2],'k-');

[X,Y,Z] = sphere(20);
X = InitR*X + InitX*ones(size(X));
Y = InitR*Y + InitY*ones(size(Y));
Z = InitR*Z + InitZ*ones(size(Z));
% surf(X,Y,Z) 
mesh(X,Y,Z)


for ii=1:NbOfLinks
%     pause;
    plot3([EstablishedLinks(ii,1)*dx, EstablishedLinks(ii,4)*dx]*1e-3,[EstablishedLinks(ii,2)*dy, EstablishedLinks(ii,5)*dy]*1e-3,[EstablishedLinks(ii,3)*dz+z_gnd, EstablishedLinks(ii,6)*dz+z_gnd]*1e-3,'b-')
    
    axis([Lx*0/4 Lx*4/4 Ly*0/4 Ly*4/4 0 Lz+z_gnd])
%     pause
end
grid
view(3)
Xp = [Lx 0 0 Lx]; 
Yp = [Ly Ly 0 0]; 
Zp = [z_gnd z_gnd z_gnd z_gnd];
patch(Xp, Yp, Zp, z_gnd,'FaceColor',gnd_color);
hold off

%-------------------------------------------------------------------------%
% Fractal Dimension                                                       %
%-------------------------------------------------------------------------%
figure(2);
clf

% Define which boundary will be reached first                             %
% NbSpheres defines the number of sphere used to draw the plot            %

NbSpheres = min([Initi ; Nx-Initi ; Initj ; Ny-Initj ; Initk ; Nz - Initk]);
fprintf('Number of spheres = %d\n', NbSpheres);
for n = 1:6
    if     (NbSpheres == Initk || NbSpheres == Nz-Initk)
        ds = dz;
    elseif (NbSpheres == Initj || NbSpheres == Ny-Initj)
        ds = dy;
    elseif (NbSpheres == Initi || NbSpheres == Nx-Initi)
        ds = dx;
    end
end

R(NbSpheres) = 0;
N(NbSpheres) = 0;
for n = 1:NbSpheres
    R(n) = InitR + n*ds;
    
    N(n) = 0;
    for m=1:NbOfLinks
        if ( (EstablishedLinks(m,4)*dx-InitX)^2 +...
             (EstablishedLinks(m,5)*dy-InitY)^2 +...
             (EstablishedLinks(m,6)*dz-InitZ)^2 <= R(n)^2 )
         N(n) = N(n)+1;
        end
    end
end

P=polyfit(log(R(3:end)), log(N(3:end)), 1);
% P=polyfit(log(R(1:3)), log(N(1:3)), 1);
plot(log(R), log(N), 'o', log(R), polyval(P, log(R)), '--')
grid
xlabel('ln(R)')
ylabel('ln(N)')
title(['Sprite fractal dimension D = ',num2str(P(1,1))])
D=P(1,1);
fprintf('Fractal dimension D = %f\n',D);
Lowest_Altitude=min(EstablishedLinks(:,6)*dz)+z_gnd;
fprintf('Lowest altitude = %f\n', Lowest_Altitude);
cd ../viz
end
%-------------------------------------------------------------------------%