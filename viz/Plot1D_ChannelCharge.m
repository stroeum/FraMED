close all
clearvars -except sims


if ~exist('sims','var') || ~isfield(sims,'pathPNGs') 
    sims = specifySimDetails();
end 

cd ../results
Q = load('CarriedCharge.dat');
if isempty(Q) || (size(Q,1)-1) == 0
    fprintf('\n*** Plot1D_ChannelCharge.m cannot be executed with current CarriedCharge.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_ChannelCharge.m script. ***\n');
end
step = (0:size(Q,1)-1)';

set(gcf,'Units','inches','OuterPosition', [20 20 20 20]/6)

plot(step,Q,'LineWidth',1,'LineStyle','-');
xlabel('step','FontSize',12);
ylabel('Q_{cha} (C)','FontSize',12);
set(gca,'FontSize',10);
axis([0 max(step) min(0) max(Q)]) ;
box on
grid on
exportgraphics(gcf,strcat(sims.pathPNGs,'/ChargeTransfer_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);

cd ../viz