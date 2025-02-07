close all
clearvars -except sims
beep  off

if ~exist('sims','var') || ~isfield(sims,'pathPNGs')
    sims = specifySimDetails();
end 
fprintf('\n*** Executing Plot2D_Fields3x3.m script. ***\n');

cd ../results/
%% Load
EthPositive = load('EthPositive.dat')*1e-5;
EthNegative = load('EthNegative.dat')*1e-5;
EzBF        = load('EnumBF.dat')*1e-5;
ExAF3D      = load('Ex3d.dat')*1e-5;
EyAF3D      = load('Ey3d.dat')*1e-5;
EzAF3D      = load('Ez3d.dat')*1e-5;
dxyz        = load('dxyz.dat')*1e-3;
Nxyz        = load('Nxyz.dat');
z_gnd       = load('z_gnd.dat')*1e-3;

%% Derive main sims
Nx = Nxyz(1);
Ny = Nxyz(2);
Nz = Nxyz(3);

dx = dxyz(1);           % _m
dy = dxyz(2);           % _m
dz = dxyz(3);           % _m

Lx = (Nx-1)*dx;         % _m
Ly = (Ny-1)*dy;         % _m
Lz = (Nz-1)*dz;         % _m

x = (0:Nx-1)'*dx;       % _m
y = (0:Ny-1)'*dy;       % _m
z = (0:Nz-1)'*dz+z_gnd; % _m
clear dxyz
cd ../viz/

%% Convert
ExAF3D  = ConvertTo3d(ExAF3D,Nxyz);
EyAF3D  = ConvertTo3d(EyAF3D,Nxyz);
EzAF3D  = ConvertTo3d(EzAF3D,Nxyz);
clear Nxyz

EzAF(Nz)= 0;
ExAF(Nz)= 0;
EyAF(Nz)= 0;

%% Plot
close all


ii = [floor((Nx+1)*1/4) (Nx+1)*2/4 ceil((Nx+1)*3/4)];
jj = [floor((Ny+1)*1/4) (Ny+1)*2/4 ceil((Ny+1)*3/4)];

pp = 1;
for nn=1:3
    for mm=1:3
        if pp==1 || pp==2 || pp==3
            pp = pp+6;
            subplot(3,3,pp)
            pp = pp-6;
        elseif pp==7 || pp==8 || pp==9
            pp = pp-6;
            subplot(3,3,pp)
            pp = pp+6;
        else
            subplot(3,3,pp)
        end
        for kk=1:Nz
            ExAF(kk) = ExAF3D(ii(mm),jj(nn),kk);
            EyAF(kk) = EyAF3D(ii(mm),jj(nn),kk);
            EzAF(kk) = EzAF3D(ii(mm),jj(nn),kk);
        end
        plot(EthPositive,z,'g--',EthNegative,z,'g--',EzBF,z,'b', EzAF,z,'r')
        axis([2*min(EthNegative) 2*max(EthPositive) min(z) max(z)])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        if pp<=3
            xlabel('Ez (kV/cm)','FontSize',12);
        end
        if pp==1 || pp==4 || pp==7
            ylabel('z (km)','FontSize',12);
        end
        title(['x = ',num2str(x(ii(mm))),' km; y = ',num2str(y(jj(nn))),' km'],'Fontsize',12)

        pp = pp+1;
    end
end
exportgraphics(gcf,strcat(sims.pathPNGs,'/FieldSpaceEvolution_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
clear mm nn pp ii jj kk