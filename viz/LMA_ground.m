function []=LMA_ground(h,panels,sims)
global alt_hist FocusArea Ground
% LMA ground plane %
if Ground.Line.Width~=0
    if strcmp(panels,'yz') || strcmp(panels,'all')
        subplot(h.yz.fig);
        hold on; plot([sims.domain.gnd, sims.domain.gnd], [FocusArea.y(1) FocusArea.y(2)], 'k', 'LineWidth', Ground.Line.Width);
    end
    if strcmp(panels,'xz') || strcmp(panels,'all')
        subplot(h.xz.fig);
        hold on; plot([FocusArea.x(1) FocusArea.x(2)], [sims.domain.gnd, sims.domain.gnd], 'k', 'LineWidth', Ground.Line.Width);
    end
    if strcmp(panels,'nz') || strcmp(panels,'all')
        subplot(h.nz.fig);
        hold on; plot([0 alt_hist.max], [sims.domain.gnd, sims.domain.gnd], 'k', 'LineWidth', Ground.Line.Width);
    end
end
