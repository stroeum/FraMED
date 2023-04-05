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
Es = load('EsEnergy.dat')*1e-9;
if isempty(Es)
    fprintf('\n*** Plot1D_EsEnergy.m cannot be executed with current EsEnergy.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_EsEnergy.m script. ***\n');
end
DischargeDipoleMoment
step = (0:size(Es)-1)';

set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

plot(step,Es,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('\epsilon_{es} (GJ)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(Es)]) ;
box on
grid on
exportgraphics(gcf,[sims.pathPNGs,'/EsEnergy_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);


cd ../viz/