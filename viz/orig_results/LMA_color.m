% LMA link color scheme
if strcmp(Links.Scheme,'Colormap')
    Links.Color      = colormap(jet(Nb.Points));
elseif strcmp(Links.Scheme,'B/W')
    Links.Color      = [0 0 0];
elseif strcmp(Links.Scheme,'R/B')
    Links.Color(1,:) = [0 0 0];
    Links.Color(2,:) = [1 0 0];
    Links.Color(3,:) = [0 0 1];
elseif strcmp(Links.Scheme,'B/R')
    Links.Color(1,:) = [0 0 0];
    Links.Color(2,:) = [0 0 1];
    Links.Color(3,:) = [1 0 0];
end
close all