clearvars -except sims
close all

% For formatting exported image size:
positionWidth = 600;
positionHeight = 800;

%% Compile and Run C++ code
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
    prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
    sims.objectType = input(prompt2,'s');
    while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
        fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
        sims.objectType = input(prompt2,'s');
    end

    % Settings to ensure proper directory referencing:
    sims.pathPNGs = ['../Figures/',sims.objectName,'/',sims.objectType,'/PNGs'];
    if ~exist(sims.pathPNGs,'dir')
        mkdir(sims.pathPNGs);
    end
    sims.pathVideos = ['../Figures/',sims.objectName,'/',sims.objectType,'/Videos'];
    if ~exist(sims.pathVideos,'dir')
        mkdir(sims.pathVideos);
    end
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

dz = dxyz(3);
Nz = Nxyz(3);
Lz = (Nz-1)*dz;
z  = (0:Nz-1)*dz+z_gnd;
clear dxyz Nxyz
if length(EthNegative)~=length(EzNumBF)
    fprintf('\n*** Plot1D_Fields.m cannot be executed with current EthNegative.dat and EnumBF.dat files.\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_Fields.m script. ***\n');
end

% Initiation requirements plot:
figure;
subplot(121)
hold on
plot(...
    EzNumBF(:),z*1e-3,'r-.',...
    1/1000*phiNumBF(:), z*1e-3, 'b-',...
    EthNegative(:),z*1e-3, 'g--',...
    EthPositive(:),z*1e-3, 'g--'...
    );
plot([min(EthNegative) max(EthPositive)],[z_gnd,z_gnd]*1e-3,'k');
% semilogx(...
%     abs(EzNum(:)), [0:size(EzNum)-1]*dz,'r-',...
%     abs(Einitiation(:)),[0:size(Einitiation)-1]*dz, 'b--',...
%     abs(EthNegative(:)),[0:size(EthNegative)-1]*dz, 'g--',...
%     abs(EthPositive(:)),[0:size(EthPositive)-1]*dz, 'g--'...
%     );

axis([min(EthNegative) max(EthPositive) 0 (Lz+z_gnd)*1e-3]);
legend('E_z [BF]','\phi [BF]','E_{th}^\pm');
xlabel('E_z (V/m) & \phi (kV)');
ylabel('Altitude (km)');
title('E_z and \phi  before the flash [BF]');
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
box on
hold off

subplot(122)
hold on
plot(...
    EzNumAF(:)*10^-9,z*1e-3,'r-.',...
    1/1000*phiNumAF(:), z*1e-3, 'b-',...
    EthNegative(:),z*1e-3, 'g--',...
    EthPositive(:),z*1e-3, 'g--'...
    );
plot([min(EthNegative) max(EthPositive)],[z_gnd,z_gnd]*1e-3,'k');
% semilogx(...
%     abs(EzNum(:)), [0:size(EzNum)-1]*dz,'r-',...
%     abs(Einitiation(:)),[0:size(Einitiation)-1]*dz, 'b--',...
%     abs(EthNegative(:)),[0:size(EthNegative)-1]*dz, 'g--',...
%     abs(EthPositive(:)),[0:size(EthPositive)-1]*dz, 'g--'...
%     );

axis([min(EthNegative) max(EthPositive) 0 (Lz+z_gnd)*1e-3]);
legend('E_z [BF]','\phi [BF]','E_{th}^\pm');
xlabel('E_z (V/m) & \phi (kV)');
ylabel('Altitude (km)','FontSize',16);
title('E_z and \phi  after the flash [AF]');
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
box on
hold off
exportgraphics(gcf,[sims.pathPNGs,'/InitiationRequirements_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);

% Potential before and after flash plot:
linewidth = 1;
figure;
%set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)
hold on
plot(phiNumBF(:)*1e-6, z*1e-3, 'b-','LineWidth',linewidth)
plot(phiNumAF(:)*1e-6, z*1e-3, 'b--','LineWidth',linewidth)
plot([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6)],[z_gnd,z_gnd]*1e-3,'k','LineWidth',linewidth);
hold off
axis([1.1*min(phiNumBF(:)*1e-6) 1.1*max(phiNumBF(:)*1e-6) 0 (Lz+z_gnd)*1e-3]);
xlabel('$\phi$ (MV)','FontSize',18,'Interpreter','latex');
ylabel('Altitude (km)','FontSize',18,'Interpreter','latex');
set(gca,'FontSize',14);
set(gca,'XMinorTick','on','YMinorTick','on')
box on
% title('E_z and \phi  before the flash [BF] and after [AF]');
grid on;
hold off
exportgraphics(gcf,[sims.pathPNGs,'/phi_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);

% Electric field thresholds plot:
fig = figure;
%set(gcf,'Position',[0,0,600,1000]);
set(gcf,'Position',[0,0,positionWidth,positionHeight]);
set(gcf,'Resize','off')
hold on
plot(EzNumBF(:)*1e-5, z*1e-3, 'r-','LineWidth',linewidth,'DisplayName','Electric Field')
plot(EthNegative(:)*1e-5,z*1e-3, 'g-.','LineWidth',linewidth,'HandleVisibility','off')
plot(EthPositive(:)*1e-5,z*1e-3, 'g-.','LineWidth',linewidth,'DisplayName','Propagation Threshold')
plot(Einitiation(:)*1e-5,z*1e-3, 'b--','LineWidth',linewidth,'DisplayName','Initiation Threshold')
plot(-Einitiation(:)*1e-5,z*1e-3, 'b--','LineWidth',linewidth,'HandleVisibility','off')
hold off
axis([-1.25*max(Einitiation*1e-5) 1.25*max(Einitiation*1e-5) 0 (Lz+z_gnd)*1e-3]);
legend('Interpreter','latex','FontSize',16,'location','northeast')
set(gca,'FontSize',16);
set(gca,'XMinorTick','on','YMinorTick','on','Tickdir','out','TickLabelInterpreter','latex')
box on
title('Electric Field Thresholds','Interpreter','latex','FontSize',28);
xlabel('$E_z$ (kV/cm)','FontSize',20,'Interpreter','latex');
ylabel('Altitude (km)','FontSize',20,'Interpreter','latex');
grid on
set(fig,'Position',[0,0,positionWidth,positionHeight]);
exportgraphics(gcf,[sims.pathPNGs,'/E_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);
%xlim([-1000000,1000000])
fprintf('\tAverage field reduction: %f\n',mean(abs(EzNumAF)./abs(EzNumBF)));

cd ../viz