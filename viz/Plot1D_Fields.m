clearvars -except sims
close all

%% Compile and Run C++ code
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    specifySimDetails;
end 
cd ../results

%% Plot results
dxyz        = load('dxyz.dat');
Nxyz        = load('Nxyz.dat');
EzNumBF     = load('EnumBF.dat');
phiNumBF    = load('phiNumBF.dat');
EzNumAF     = load('EnumAF.dat');
phiNumAF    = load('phiNumAF.dat');
Einitiation = load('Einitiation.dat');
EthPositive = load('EthPositive.dat');
EthNegative = load('EthNegative.dat');
z_gnd       = load('z_gnd.dat');

if length(EthNegative)~=length(EzNumBF)
    fprintf('\n*** Plot1D_Fields.m cannot be executed with current EthNegative.dat and EnumBF.dat files.\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_Fields.m script. ***\n');
end
cd ../viz

dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz

% Determine the magnitude of the values:
spatialFactor = checkMagnitude(z(:));
fieldFactor = checkMagnitude([EzNumBF(:); EzNumAF(:); Einitiation(:); EthPositive(:); EthNegative(:)]);
potentialFactor = checkMagnitude([phiNumBF(:); phiNumAF(:)]);
maxfieldLim = 1.1*max([fieldFactor.Number*max([-min(EzNumBF) max(EzNumBF) -min(EzNumAF) max(EzNumAF) -min(EthNegative) max(EthPositive) max(Einitiation)]) potentialFactor.Number*max([-min(phiNumBF) max(phiNumBF) -min(phiNumAF) max(phiNumAF)])]);

% Initiation requirements plot:
figure;
subplot(121)
hold on
plot(...
    fieldFactor.Number*EzNumBF(:),z*spatialFactor.Number,'r-.',...
    potentialFactor.Number*phiNumBF(:), z*spatialFactor.Number, 'b-',...
    fieldFactor.Number*EthNegative(:),z*spatialFactor.Number, 'g--',...
    fieldFactor.Number*EthPositive(:),z*spatialFactor.Number, 'g--'...
    );
plot([fieldFactor.Number*min(EthNegative) fieldFactor.Number*max(EthPositive)],[z_gnd,z_gnd]*spatialFactor.Number,'k');

axis([-maxfieldLim maxfieldLim 0 (Lz+z_gnd)*spatialFactor.Number]);
legend('$E_z$ [AF]','$\phi$ [AF]','$E_{th}^\pm$','interpreter','latex');
xlabel(strcat('$E_z$ (',fieldFactor.Unit,'V/m) \& $\phi$ (',potentialFactor.Unit,'V)'),'interpreter','latex');
ylabel(strcat('Altitude (',spatialFactor.Unit,'m)'),'interpreter','latex');
title('$E_z$ and $\phi$ before the flash [BF]','interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')
grid on
box on
hold off

subplot(122)
hold on
plot(...
    EzNumAF(:)*fieldFactor.Number,z*spatialFactor.Number,'r-.',...
    potentialFactor.Number*phiNumAF(:), z*spatialFactor.Number, 'b-',...
    fieldFactor.Number*EthNegative(:),z*spatialFactor.Number, 'g--',...
    fieldFactor.Number*EthPositive(:),z*spatialFactor.Number, 'g--'...
    );
plot([fieldFactor.Number*min(EthNegative) fieldFactor.Number*max(EthPositive)],[z_gnd,z_gnd]*spatialFactor.Number,'k');
axis([-maxfieldLim maxfieldLim 0 (Lz+z_gnd)*spatialFactor.Number]);
legend('$E_z$ [AF]','$\phi$ [AF]','$E_{th}^\pm$','interpreter','latex');
xlabel(strcat('$E_z$ (',fieldFactor.Unit,'V/m) \& $\phi$ (',potentialFactor.LaTeX,'V)'),'interpreter','latex');
ylabel(strcat('Altitude (',spatialFactor.LaTeX,'m)'),'interpreter','latex');
title('$E_z$ and $\phi$ after the flash [AF]','interpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')
grid on
box on
hold off
exportgraphics(gcf,strcat(sims.pathPNGs,'/InitiationRequirements_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);

% Potential before and after flash plot:
linewidth = 1;
figure;
%set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(phiNumBF(:)*potentialFactor.Number, z*spatialFactor.Number, 'b-','LineWidth',linewidth)
plot(phiNumAF(:)*potentialFactor.Number, z*spatialFactor.Number, 'b--','LineWidth',linewidth)
plot([1.1*min(phiNumBF(:)*potentialFactor.Number) 1.1*max(phiNumBF(:)*potentialFactor.Number)],[z_gnd,z_gnd]*spatialFactor.Number,'k','LineWidth',linewidth);
hold off
axis([1.1*min(phiNumBF(:)*potentialFactor.Number) 1.1*max(phiNumBF(:)*potentialFactor.Number) 0 (Lz+z_gnd)*spatialFactor.Number]);
legend('$\phi$ [BF]','$\phi$ [AF]','interpreter','latex','location','best');
xlabel(strcat('$\phi$ (',potentialFactor.LaTeX,'V)'),'FontSize',18,'Interpreter','latex');
ylabel(strcat('Altitude (',spatialFactor.LaTeX,'m)'),'FontSize',18,'Interpreter','latex');
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on','TickLabelInterpreter','latex')
box on
title('$\phi$  before the flash [BF] and after [AF]','interpreter','latex');
grid on;
hold off
exportgraphics(gcf,strcat(sims.pathPNGs,'/phi_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);

% Electric field thresholds plot:
fig = figure;
set(gcf,'Position',[0,0,600,800]);
set(gcf,'Resize','off')
hold on
plot(EzNumBF(:)*fieldFactor.Number, z*spatialFactor.Number, 'r-','LineWidth',linewidth,'DisplayName','Electric Field')
plot(EthNegative(:)*fieldFactor.Number,z*spatialFactor.Number, 'g-.','LineWidth',linewidth,'HandleVisibility','off')
plot(EthPositive(:)*fieldFactor.Number,z*spatialFactor.Number, 'g-.','LineWidth',linewidth,'DisplayName','Propagation Threshold')
plot(Einitiation(:)*fieldFactor.Number,z*spatialFactor.Number, 'b--','LineWidth',linewidth,'DisplayName','Initiation Threshold')
plot(-Einitiation(:)*fieldFactor.Number,z*spatialFactor.Number, 'b--','LineWidth',linewidth,'HandleVisibility','off')
hold off
axis([1.25*fieldFactor.Number*min([-max(Einitiation) -max(EzNumBF) min(EzNumBF)]) 1.25*fieldFactor.Number*max([max(Einitiation) -min(EzNumBF) max(EzNumBF)]) 0 (Lz+z_gnd)*spatialFactor.Number]);
legend('Interpreter','latex','FontSize',16,'location','northeast')
set(gca,'FontSize',16,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
box on
title('Electric Field Thresholds','Interpreter','latex','FontSize',28);
xlabel(strcat('$E_z$ (',fieldFactor.LaTeX,'V/m)'),'FontSize',20,'Interpreter','latex');
ylabel(strcat('Altitude (',spatialFactor.LaTeX,'m)'),'FontSize',20,'Interpreter','latex');
grid on
set(gcf,'Position',[0,0,600,800]);
exportgraphics(gcf,strcat(sims.pathPNGs,'/E_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
fprintf('\tAverage field reduction: %f\n',mean(abs(EzNumAF)./abs(EzNumBF)));