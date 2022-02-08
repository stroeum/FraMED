function []=LMA_labels(h,panels)
global alt_hist FocusArea Font Nb
% LMA axes %

if strcmp(panels,'xy') || strcmp(panels,'all')
    subplot(h.xy.fig);
    h.xy.label.x.textbox     = xlabel(h.xy.label.x.string,'Units','Points');
    h.xy.label.x.position    = get(h.xy.label.x.textbox,'Position');
    h.xy.label.x.position(2) = -h.xy.label.spacing.x;
    set(h.xy.label.x.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.xy.label.x.position)
    h.xy.label.y.textbox     = ylabel(h.xy.label.y.string,'Units','Points');
    h.xy.label.y.position    = get(h.xy.label.y.textbox,'Position');
    h.xy.label.y.position(1) = -h.xy.label.spacing.y;
    set(h.xy.label.y.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.xy.label.y.position)

end

if strcmp(panels,'yz') || strcmp(panels,'all')
    subplot(h.yz.fig);
    h.yz.label.x.textbox     = xlabel(h.yz.label.x.string,'Units','Points');
    h.yz.label.x.position    = get(h.yz.label.x.textbox,'Position');
    h.yz.label.x.position(2) = -h.yz.label.spacing.x;
    set(h.yz.label.x.textbox, 'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.yz.label.x.position)
    h.yz.label.y.textbox     = ylabel(h.yz.label.y.string,'Units','Points');
    h.yz.label.y.position    = get(h.yz.label.y.textbox,'Position');
    h.yz.label.y.position(1) = -h.yz.label.spacing.y;
    set(h.yz.label.y.textbox, 'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.yz.label.y.position)
end

if strcmp(panels,'xz') || strcmp(panels,'all')
    subplot(h.xz.fig);
    h.xz.label.x.textbox     = xlabel(h.xz.label.x.string,'Units','Points');
    h.xz.label.x.position    = get(h.xz.label.x.textbox,'Position');
    h.xz.label.x.position(2) = -h.xz.label.spacing.x;
    set(h.xz.label.x.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.xz.label.x.position)
    h.xz.label.y.textbox     = ylabel(h.xz.label.y.string,'Units','Points');
    h.xz.label.y.position    = get(h.xz.label.y.textbox,'Position');
    h.xz.label.y.position(1) = -h.xz.label.spacing.y;
    set(h.xz.label.y.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.xz.label.y.position)
end

if strcmp(panels,'nz') || strcmp(panels,'all')
    subplot(h.nz.fig);
    h.nz.label.x.textbox     = xlabel(h.nz.label.x.string,'Units','Points');
    h.nz.label.x.position    = get(h.nz.label.x.textbox,'Position');
    h.nz.label.x.position(2) = -h.nz.label.spacing.x;
    set(h.nz.label.x.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.nz.label.x.position)
    h.nz.label.y.textbox     = ylabel(h.nz.label.y.string,'Units','Points');
    h.nz.label.y.position    = get(h.nz.label.y.textbox,'Position');
    h.nz.label.y.position(1) = -h.nz.label.spacing.y;
    set(h.nz.label.y.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.nz.label.y.position)
    h.nz.title.textbox       = title(h.nz.title.string,'Units','normalized');
    h.nz.title.position      = get(h.nz.title.textbox,'Position');
    h.nz.title.position(2)   = 1;
    set(h.nz.title.textbox,'Units','normalized','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Bottom','Position',h.nz.title.position)
end

if strcmp(panels,'yz') || strcmp(panels,'all')
    subplot(h.tz.fig);
    h.tz.label.x.textbox     = xlabel(h.tz.label.x.string,'Units','Points');
    h.tz.label.x.position    = get(h.tz.label.x.textbox,'Position');
    h.tz.label.x.position(2) = 0;
    set(h.tz.label.x.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Top','Position',h.tz.label.x.position)
    h.tz.label.y.textbox     = ylabel(h.tz.label.y.string,'Units','Points');
    h.tz.label.y.position    = get(h.tz.label.y.textbox,'Position');
    h.tz.label.y.position(1) = -h.tz.label.spacing.y;
    set(h.tz.label.y.textbox,'Units','Points','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Middle','Position',h.tz.label.y.position)
    h.tz.title.textbox       = title(h.tz.title.string,'Units','normalized');
    h.tz.title.position      = get(h.tz.title.textbox,'Position');
    h.tz.title.position(2)   = 1;
    set(h.tz.title.textbox,'Units','normalized','FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Bottom','Position',h.tz.title.position)
end

end