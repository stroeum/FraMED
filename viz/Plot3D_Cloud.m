function []=Plot3D_Cloud()
beep  off
close all
clearvars -except sims

% figure(1);
%-------------------------------------------------------------------------%
% Load data files                                                         %
%-------------------------------------------------------------------------%
cd ../results
dxyz       = load('dxyz.dat');
Nxyz       = load('Nxyz.dat');
z_gnd      = load('z_gnd.dat');
SavingStep = load('rho3dStep.dat');
Links      = load('EstablishedLinks.dat');

%-------------------------------------------------------------------------%
% Derive main parameters                                                  %
%-------------------------------------------------------------------------%
Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

X = (0:Nx-1)'*dx;       % _m
Y = (0:Ny-1)'*dy;       % _m
Z = (0:Nz-1)'*dz+z_gnd; % _m

NbOfLinks = size(Links);
NbOfLinks = NbOfLinks(1);
K = floor(NbOfLinks/SavingStep);
clear dxyz Links

fprintf('Select your vizualisation:\n');
fprintf('1. Isosurfaces\n');
fprintf('2. Slices\n');
fprintf('3. Scatter3D\n');
fprintf('exit else\n');
% choice =4;
% fprintf('choice >> %d\n',choice);
choice = 2;%input('choice >> ');

if (choice == 1)
    %---------------------------------------------------------------------%
    % Isosurface                                                          %
    %---------------------------------------------------------------------%
    rho   = load('rho3d0.dat');
    rho3D   = ConvertTo3d(rho,Nxyz);
    %     save rho3D.mat rho3D -mat
    clear Nxyz
    % load rho3D.mat -mat
    
    figure;
    % 1. Prepare the Data %
    data = rho3D;
    % data = smooth3(data,'gaussian',1);
    data = smooth3(data,'box',5);
    
    % 2. Create the Isosurface and Set Properties %
    isoval = .5;
    h = patch(isosurface(X,Y,Z,data,isoval),...
        'FaceColor','none',...
        'EdgeColor','blue',...
        'AmbientStrength',.2,...
        'SpecularStrength',.7,...
        'DiffuseStrength',.4);
    isonormals(data,h)
    
    % 3. Create the Isocaps and Set Properties %
    patch(isocaps(X,Y,Z,data,isoval),...
        'FaceColor','interp',...
        'EdgeColor','none')
    colormap hsv
    
    % 4. Define the View %
    daspect([1,1,1])
    axis([0 Lx 0 Ly 0 Lz+z_gnd])% tight
    view(3)
    
    % 5. Add Lighting %
    camlight right
    camlight left
    set(gcf,'Renderer','zbuffer');
    lighting phong
elseif(choice == 2)
    %---------------------------------------------------------------------%
    % Slice                                                               %
    %---------------------------------------------------------------------%
    rho   = load('rho3d0.dat');
    rho3D   = ConvertTo3d(rho,Nxyz);
%     save rho3D.mat rho3D -mat
    % load rho3D.mat -mat
    clear Nxyz
    
    figure;
    Record   = 0; %input('Record the movie? (1: yes, else: no)\n>> ');
    N        = Nz;
    Movie(N) = getframe;
    
    [X,Y,Z]  = meshgrid(X,Y,Z);
    size(X)
    size(Y)
    size(Z)
    size(rho3D)
    
    slice(X,Y,Z,rho3D,Lx,Ly,Lz);
    set(gca,'nextplot','replacechildren');
    axis equal tight
    if (Record ~= 1)
        for n = 1:N
            Sx      = Lx/2;
            Sy      = Ly/2;
            Sz      = [0 (n-1)*dz]+z_gnd;
            h=slice(X*1e-3,Y*1e-3,Z*1e-3,rho3D,Sx*1e-3,Sy*1e-3,Sz*1e-3,'makima');
            set(h,'edgecolor','none');
            colorbar
            xlabel('x-axis');
            ylabel('y-axis');
            zlabel('z-axis');
            title('Charge density \rho [nC]')
            axis([0 Lx 0 Ly 0 Lz+z_gnd]*1e-3)
            view([315+((n-1)*90/(N-2)) 25])
            Movie(n) = getframe;
            
        end
        colorbar;
    elseif (Record == 1)
        for n = 1:N
            Sx      = [Lx/2 Lx];
            Sy      = [Ly/2 Ly];
            Sz      = [0 (n-1)*dz]+z_gnd;
            slice(X*1e-3,Y*1e-3,Z*1e-3,rho3D,Sx*1e-3,Sy*1e-3,Sz*1e-3);
            
            colorbar;
            xlabel('x-axis');
            ylabel('y-axis');
            zlabel('z-axis');
            title('Charge density \rho [nC]')
            axis([0 Lx 0 Ly 0 Lz+z_gnd]*1e-3)
            %         view([315+((n-1)*90/(N-2)) 25]);
            view([315 25])
            Movie(n) = getframe;
        end
%         writeVide(Movie,'cloud.avi','quality',100);
    end
elseif(choice == 3)
    %---------------------------------------------------------------------%
    % Scatter3                                                            %
    %---------------------------------------------------------------------%
    rhoAmb     = load('rho3d0.dat');
    isStatic = 1;%input('Display a movie? (1: yes, else: no)\n>> ');
    
    
    if(isStatic~=1)
        rho   = rhoAmb;
        clear Nxyz
        n=1;
        for k=1:Nz
            for j=1:Ny
                for i=1:Nx
                    X(n) = (i-1)*dx;
                    Y(n) = (j-1)*dy;
                    Z(n) = (k-1)*dz+z_gnd;
                    n=n+1;
                end
            end
        end
        caxis([-15 15])
        S=2.5e3*abs(rho/max(rho))+1;
        C=rho;
        scatter3(X,Y,Z,S,C,'.')
        colorbar;
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        title('Charge density \rho [nC]');
        axis([0 Lx 0 Ly 0 Lz+z_gnd])
        view([300 5])
        %         view([0 1 0])
    else
        n        = 1;
        for k=1:Nz
            for j=1:Ny
                for i=1:Nx
                    X(n) = (i-1)*dx;
                    Y(n) = (j-1)*dy;
                    Z(n) = (k-1)*dz+z_gnd;
                    n=n+1;
                end
            end
        end
        
        Movie(K+1) =getframe();
        for k=0:K
            fname = ['rho3d',int2str(k*SavingStep),'.dat'];
            rho   = load(fname)-rhoAmb;
            S=1e3*abs(rho/max(rho))+1;
            C=rho;
            scatter3(X*1e-3,Y*1e-3,Z*1e-3,S,C,'.')
            caxis([-100 100])
            colorbar;
            xlabel('x-axis (km)','FontSize',12);
            ylabel('y-axis (km)','FontSize',12);
            zlabel('z-axis (km)','FontSize',12);
            title(['\rho  (nC) after ', int2str(k*SavingStep) ,' step(s)'],'FontSize',12);
            axis([0 Lx 0 Ly 0 Lz+z_gnd]*1e-3);
            set(gca,'FontSize',10);
            view([300 5])
            % view([0 1 0])
            axis xy equal
            Movie(k+1) = getframe(gcf,[0,0, 560, 420]);
        end
        %         movie2avi(Movie,'cloud.avi','quality',100);
    end
end
cd ../viz

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [AA] = ConvertTo3d(A,B)
    [M, N] = size(A);
    AA = zeros(B');
    for m=1:M
        for n=1:N
            ii = rem(m,B(1));
            if(ii==0)
                ii = B(1);
            end
            jj = n;
            kk = (m-ii)/B(1)+1;
            AA(ii,jj,kk) = A(m,n);
        end
    end
end



% function [AA] = ConvertTo1d(A,B)
% [M N] = size(A);
% p=1;
% for m=1:M
%     for n=1:N
%         ii = rem(m,B(1));
%         if(ii==0)
%             ii =B(1);
%         end
%         jj = n;
%         kk = (m-ii)/B(1)+1;
%         AA(p,1) = ii;
%         AA(p,2) = jj;
%         AA(p,3) = kk;
%         AA(p,4) = A(m,n);
%         p = p+1;
%     end
% end
%
% end
%-------------------------------------------------------------------------%