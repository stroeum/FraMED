function []= Plot3D_Tree()
close all
clear all
clc
cd ../bin/results/

%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%

load dxyz.dat               -ascii
load Nxyz.dat               -ascii
load InitPoint.dat          -ascii
load EstablishedLinks.dat   -ascii
load z_gnd.dat              -ascii

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
% Map ColorScale                                                          %
%-------------------------------------------------------------------------%
gnd_color = [.5 .5 .5];
color     = colormap(jet(NbOfLinks));
isMonochrome = 1;%input('Is plot Monochrome? (1: yes, else: no)\n>> ');
if (isMonochrome == 1)
    for ii=1:NbOfLinks
        color(ii,:) = [0 0 1];
    end
end

%-------------------------------------------------------------------------%
% Draw the tree                                                           %
%-------------------------------------------------------------------------%
figure(1);
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/3)
hold 
%Cloud3D(3e3, 18e3, Lx, Ly, Lz, z_gnd);

Xp = [Lx 0 0 Lx]*1e-3;
Yp = [Ly Ly 0 0]*1e-3;
Zp = [z_gnd z_gnd z_gnd z_gnd]*1e-3;
patch(Xp, Yp, Zp, z_gnd,'FaceColor',gnd_color);



% xc = InitX;
% yc = InitY;
% zc = InitZ;
% lx = 4*dx;%0.0025;%2e-3;
% ly = 4*dy;%0.0025;%2e-3;
% lz = 4*dz;%0.0025;%2e-3;

% [X,Y,Z] = sphere(20);
% X = InitR*X + InitX*ones(size(X));
% Y = InitR*Y + InitY*ones(size(Y));
% Z = InitR*Z + InitZ*ones(size(Z));
% % surf(X,Y,Z)
% mesh(X,Y,Z)
% Movie(NbOfLinks) = getframe(gcf,[0,0, 560, 420]);
for ii=1:NbOfLinks
    if(rem(ii,10)==0)
        %         pause;
    end
       
    plot3(...
        [EstablishedLinks(ii,1)*dx      , EstablishedLinks(ii,4)*dx      ]*1e-3,...
        [EstablishedLinks(ii,2)*dy      , EstablishedLinks(ii,5)*dy      ]*1e-3,...
        [EstablishedLinks(ii,3)*dz+z_gnd, EstablishedLinks(ii,6)*dz+z_gnd]*1e-3,...
        'Color',color(ii,:));
    axis([Lx*0/4 Lx*4/4 Ly*0/4 Ly*4/4 z_gnd 2/2*(Lz+z_gnd)]*1e-3)
%     axis equal
    xlabel('x (km)','FontSize',12);
    ylabel('y (km)','FontSize',12);
    zlabel('z (km)','FontSize',12);
    title(['Lightning discharge after ', int2str(ii) ,' step(s)'],'FontSize',12,'FontWeight','bold');
    set(gca,'FontSize',10);
    view([45 10])
%     Movie(ii) = getframe(gcf,[0,0, 599, 553]);
%     Movie(ii) = getframe;
end
grid
Movie(NbOfLinks+1) = getframe(gcf,[0,0, 599, 553]);
hold off;

Record = input('Record the movie? (1: yes, else: no)\n>> ');
if (Record == 1)
%     close all;
%     movie(Movie); % play the movie
    movie2avi(Movie,'lightning.avi','quality',100);
end
cd ../../viz

end
%-------------------------------------------------------------------------%
