% Initiation of discharge tree
% Figure.Outerposition = get
if strcmp(Figure.Output,'Movie')
    Movie(1) = getframe(h.fig);
    Movie(Nb.Points+1) = getframe(h.fig);
end

alt_hist = LMA_alt_hist(alt_hist,0);

subplot(h.tz.fig);
hold on; plot(0,(Links.data(1,3)*d.z+z_gnd),'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

subplot(h.nz.fig);
reset(subplot(h.nz.fig));
if strcmp(Links.Type,'Both')
    h.nz.note.string = [int2str(max(alt_hist.data)),' pts '];
    if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
        plot(alt_hist.data,z,'k-');
    elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        plot(alt_hist.upper,z,'Color',Links.Color(2,:));
        hold on;
        plot(alt_hist.lower,z,'Color',Links.Color(3,:));
        plot(alt_hist.data,z,'k-');
        hold off;
    end
elseif strcmp(Links.Type,'Upper')
    h.nz.note.string = [int2str(max(alt_hist.upper)),' pts '];
    if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
        plot(alt_hist.data,z,'k-');
    elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        plot(alt_hist.upper,z,'Color',Links.Color(2,:));
    end
elseif strcmp(Links.Type,'Lower')
    h.nz.note.string = [int2str(max(alt_hist.lower)),' pts '];
    if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
        plot(alt_hist.data,z,'k-');
    elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        plot(alt_hist.lower,z,'Color',Links.Color(3,:));
    end
end
h.nz.note.string = [int2str(1),' pts'];
text(alt_hist.max,FocusArea.z(2),h.nz.note.string,'FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Top','HorizontalAlignment','Right')
LMA_ground(h,'nz')
LMA_labels(h,'nz')
LMA_axes(h,'nz')

subplot(h.xz.fig);
hold on; plot(Links.data(1,1)*d.x,(Links.data(1,3)*d.z+z_gnd),'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

subplot(h.xy.fig);
hold on; plot(Links.data(1,1)*d.x,Links.data(1,2)*d.y,'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

subplot(h.yz.fig);
hold on; plot((Links.data(1,3)*d.z+z_gnd),Links.data(1,2)*d.y,'x-','LineWidth',Links.Line.Width,'MarkerSize',Links.MarkerSize,'MarkerEdgeColor',Links.Color(1,:),'MarkerFaceColor',Links.Color(1,:));

if strcmp(Figure.Output,'Movie')
    Movie(2) = getframe(h.fig);
end