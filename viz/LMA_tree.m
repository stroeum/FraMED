% LMA trees %
for ii=1:Nb.Links
%     pause
    alt_hist = LMA_alt_hist(alt_hist,ii,sims);

    % Current link type %
    if(Links.data(ii,6)*sims.domain.dz >= Init.z) % Upper point
        tmp.type = 'Upper';
    else % Lower point
        tmp.type = 'Lower';
    end

    % Current link color %
    if strcmp(Links.Scheme,'Colormap')
        tmp.color = Links.Color(ii,:);
    elseif strcmp(Links.Scheme,'B/W')
        tmp.color = 'k';
    elseif strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
        if strcmp(tmp.type,'Upper')
            tmp.color = Links.Color(2,:);
        else
            tmp.color = Links.Color(3,:);
        end
    end

    if strcmp(Links.Type,tmp.type) || strcmp(Links.Type,'Both')
        subplot(h.tz.fig);
        plot(ii,sims.spatialFactor.Number*(Links.data(ii,6)*sims.domain.dz+sims.domain.gnd),'s-','LineWidth',Links.Line.Width,'MarkerSize',floor(Links.MarkerSize/7),'MarkerEdgeColor',tmp.color,'MarkerFaceColor',tmp.color);

        subplot(h.nz.fig);
        reset(subplot(h.nz.fig))
        if strcmp(Links.Type,'Both')
            h.nz.note.string = [int2str(sum(alt_hist.data)),' pts '];
            if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
                plot(alt_hist.data,z,'k-');
            elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
                plot(alt_hist.upper,z,'Color',Links.Color(2,:));
                hold on;
                plot(alt_hist.lower,z,'Color',Links.Color(3,:));
                plot(alt_hist.data,z,'k--');
                hold off;
            end
        elseif strcmp(Links.Type,'Upper')
            h.nz.note.string = [int2str(sum(alt_hist.upper)),' pts '];
            if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
                plot(alt_hist.upper,z,'k-');
            elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
                plot(alt_hist.upper,z,'Color',Links.Color(2,:));
            end
        elseif strcmp(Links.Type,'Lower')
            h.nz.note.string = [int2str(sum(alt_hist.lower)),' pts '];
            if strcmp(Links.Scheme,'Colormap') || strcmp(Links.Scheme,'B/W')
                plot(alt_hist.lower,z,'k-');
            elseif  strcmp(Links.Scheme,'R/B') || strcmp(Links.Scheme,'B/R')
                plot(alt_hist.lower,z,'Color',Links.Color(3,:));
            end
        end
        LMA_ground(h,'nz',sims)
        LMA_labels(h,'nz')
        LMA_axes(h,'nz')
        text(alt_hist.max,FocusArea.z(2),h.nz.note.string,'FontName',Font.Name,'Fontsize',Font.Size.Labels,'VerticalAlignment','Top','HorizontalAlignment','Right')

        sf1 = subplot(h.xz.fig);
        plot(sims.spatialFactor.Number*[Links.data(ii,1)*sims.domain.dx, Links.data(ii,4)*sims.domain.dx],sims.spatialFactor.Number*[Links.data(ii,3)*sims.domain.dz+sims.domain.gnd, Links.data(ii,6)*sims.domain.dz+sims.domain.gnd],'Color',tmp.color,'LineWidth',Links.Line.Width);   
        
        sf2 = subplot(h.yz.fig);
        plot(sims.spatialFactor.Number*[(Links.data(ii,3)*sims.domain.dz)+sims.domain.gnd, (Links.data(ii,6)*sims.domain.dz)+sims.domain.gnd],sims.spatialFactor.Number*[Links.data(ii,2)*sims.domain.dy, Links.data(ii,5)*sims.domain.dy],'Color',tmp.color,'LineWidth',Links.Line.Width);

        sf3 = subplot(h.xy.fig);
        plot(sims.spatialFactor.Number*[Links.data(ii,1)*sims.domain.dx, Links.data(ii,4)*sims.domain.dx],sims.spatialFactor.Number*[Links.data(ii,2)*sims.domain.dy, Links.data(ii,5)*sims.domain.dy],'Color',tmp.color,'LineWidth',Links.Line.Width);
    end
    clear tmp

    if strcmp(Figure.Output,'Movie')
        Movie = getframe(h.fig);
        writeVideo(vFile,Movie);
    end
end

clear ii