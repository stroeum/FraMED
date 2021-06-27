function []=LMA_ground(h,panels)
global alt_hist FocusArea Ground z_gnd
% LMA ground plane %
if Ground.Line.Width~=0
    if strcmp(panels,'yz') || strcmp(panels,'all')
        subplot(h.yz.fig);
        hold on; plot([z_gnd, z_gnd], [FocusArea.y(1) FocusArea.y(2)], 'k', 'LineWidth', Ground.Line.Width);
    end
    if strcmp(panels,'xz') || strcmp(panels,'all')
        subplot(h.xz.fig);
        hold on; plot([FocusArea.x(1) FocusArea.x(2)], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);
    end
    if strcmp(panels,'nz') || strcmp(panels,'all')
        subplot(h.nz.fig);
        hold on; plot([0 alt_hist.max], [z_gnd, z_gnd], 'k', 'LineWidth', Ground.Line.Width);
    end
end
