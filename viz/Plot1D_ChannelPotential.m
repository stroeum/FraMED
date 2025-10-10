close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs')
    specifySimDetails;
end 

cd ../results/
% 1 pica = 12 points = 1/6 inch
%set(gca,'Units','inches','Position',[2 2 20 15]/6,'OuterPosition', [0 0 30 25]/6)
set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

phio = load('ChannelPotentials.dat');
if isempty(phio) || size(phio,1)-1 == 0
    fprintf('\n*** Plot1D_ChannelPotential.m cannot be executed with current ChannelPotentials.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_ChannelPotential.m script. ***\n');
end
cd ../viz/

step = (0:(size(phio,1)-1))';
potentialFactor = checkMagnitude(phio);
plot(step,potentialFactor.Number*phio,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel(strcat('\phi_0 (',potentialFactor.Unit,'V)'),'FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) potentialFactor.Number*min(phio) potentialFactor.Number*max(phio)]) ;
box on
grid on
exportgraphics(gcf,strcat(sims.pathPNGs,'/ChannelPotential_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
