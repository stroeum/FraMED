close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs')
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
end 

cd ../results/
% 1 pica = 12 points = 1/6 inch
%set(gca,'Units','inches','Position',[2 2 20 15]/6,'OuterPosition', [0 0 30 25]/6)
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

phio = load('ChannelPotentials.dat')*1e-6;
if isempty(phio)
    fprintf('\n*** Plot1D_ChannelPotential.m cannot be executed with current ChannelPotentials.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_ChannelPotential.m script. ***\n');
end
step = (0:size(phio)-1)';

plot(step,phio,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('\phi_0 (MV)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(phio) max(phio)]) ;
box on
grid on
exportgraphics(gcf,[sims.pathPNGs,'/ChannelPotential_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);

cd ../viz/