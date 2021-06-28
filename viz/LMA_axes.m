function []=LMA_axes(h,panels)
global alt_hist FocusArea Font Nb
% LMA axes %

if strcmp(panels,'xy') || strcmp(panels,'all')
    subplot(h.xy.fig);
    set(h.xy.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name,'XMinorTick','on','YMinorTick','on')
    axis([FocusArea.x(1) FocusArea.x(2) FocusArea.y(1) FocusArea.y(2)]);
    box on
end

if strcmp(panels,'yz') || strcmp(panels,'all')
    subplot(h.yz.fig);
    set(h.yz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name,'XMinorTick','on','YMinorTick','on')
    % set(h.yz.fig,'YTickLabel','');
    axis([FocusArea.z(1) FocusArea.z(2) FocusArea.y(1) FocusArea.y(2)]);
    box on
end

if strcmp(panels,'xz') || strcmp(panels,'all')
    subplot(h.xz.fig);
    set(h.xz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name,'XMinorTick','on','YMinorTick','on')
    % set(h.xz.fig,'XTickLabel','');
    axis([FocusArea.x(1) FocusArea.x(2) FocusArea.z(1) FocusArea.z(2)]);
    box on
end

if strcmp(panels,'nz') || strcmp(panels,'all')
    subplot(h.nz.fig);
    set(h.nz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name,'XMinorTick','on','YMinorTick','on')
    axis([0 alt_hist.max FocusArea.z(1) FocusArea.z(2)]);
    box on
end

if strcmp(panels,'yz') || strcmp(panels,'all')
    subplot(h.tz.fig);
    set(h.tz.fig,'FontSize',Font.Size.Axis,'FontName',Font.Name,'XMinorTick','on','YMinorTick','on')
    axis([0 Nb.Points FocusArea.z(1) FocusArea.z(2)]);
    box on
end

end