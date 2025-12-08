close all
clearvars -except sims
beep  off

if ~exist('sims','var')
    specifySimDetails;
end 
fprintf('\n*** Executing Plot2D_Fields3x3.m script. ***\n');

cd ../results/
%% Load
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
EzBF        = load('EnumBF.dat');
ExAF3D      = load('Ex3d.dat');
EyAF3D      = load('Ey3d.dat');
EzAF3D      = load('Ez3d.dat');

%% Derive main sims
x        = (0:sims.domain.dx:sims.domain.maxx)';
y        = (0:sims.domain.dy:sims.domain.maxy)';
z        = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)';
cd ../viz/

% Determine the magnitude of the values:
fieldFactor = checkMagnitude([EzBF(:); ExAF3D(:); EyAF3D(:); EzAF3D(:); EthPositive(:); EthNegative(:)]);
maxfieldLim = 1.1*fieldFactor.Number*max([max(max(abs(EzAF3D))) -min(EzBF) max(EzBF) -min(EthNegative) max(EthPositive)]);

%% Convert
ExAF3D  = convertTo3d(ExAF3D,sims);
EyAF3D  = convertTo3d(EyAF3D,sims);
EzAF3D  = convertTo3d(EzAF3D,sims);

EzAF(sims.domain.Nz)= 0;
ExAF(sims.domain.Nz)= 0;
EyAF(sims.domain.Nz)= 0;

%% Plot
close all


ii = [floor((sims.domain.Nx+1)*1/4) (sims.domain.Nx+1)*2/4 ceil((sims.domain.Nx+1)*3/4)];
jj = [floor((sims.domain.Ny+1)*1/4) (sims.domain.Ny+1)*2/4 ceil((sims.domain.Ny+1)*3/4)];

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
        for kk=1:sims.domain.Nz
            ExAF(kk) = ExAF3D(ii(mm),jj(nn),kk);
            EyAF(kk) = EyAF3D(ii(mm),jj(nn),kk);
            EzAF(kk) = EzAF3D(ii(mm),jj(nn),kk);
        end
        plot(fieldFactor.Number*EthPositive,sims.spatialFactor.Number*z,'g--',fieldFactor.Number*EthNegative,sims.spatialFactor.Number*z,'g--',fieldFactor.Number*EzBF,sims.spatialFactor.Number*z,'b', fieldFactor.Number*EzAF,sims.spatialFactor.Number*z,'r')
        axis([-maxfieldLim maxfieldLim sims.spatialFactor.Number*min(z) sims.spatialFactor.Number*max(z)])
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        if pp<=3
            xlabel(strcat('Ez (',fieldFactor.Unit,'V/m)'),'FontSize',12);
        end
        if pp==1 || pp==4 || pp==7
            ylabel(strcat('z (',sims.spatialFactor.Unit,'m)'),'FontSize',12);
        end
        title(strcat("x = ",num2str(sims.spatialFactor.Number*x(ii(mm)))," ",sims.spatialFactor.Unit,"m; y = ",num2str(sims.spatialFactor.Number*y(jj(nn)))," ",sims.spatialFactor.Unit,'m'),'Fontsize',12)

        pp = pp+1;
    end
end
exportgraphics(gcf,strcat(sims.pathPNGs,'/FieldSpaceEvolution_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
clear mm nn pp ii jj kk