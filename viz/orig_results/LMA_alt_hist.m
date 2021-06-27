function [alt_hist]=LMA_alt_hist(alt_hist,step)
% Altitude histogram %
global d Init Links N Nb

if step == -1
    alt_hist.tmp(N.z)    = 0;
    alt_hist.tmp(Links.data(1,3)+1) = 1;
    for ii = 1:Nb.Links
        alt_hist.tmp(Links.data(ii,6)+1) = alt_hist.tmp(Links.data(ii,6)+1) +1;
    end
    alt_hist.max = max(alt_hist.tmp);
elseif step == 0
    alt_hist.data(1:N.z)              = 0;
    alt_hist.upper(1:N.z)             = 0;
    alt_hist.lower(1:N.z)             = 0;
    alt_hist.data(Links.data(1,3)+1)  = 1;
    if(Links.data(1,3)*d.z >= Init.z)
        alt_hist.upper(Links.data(1,3)+1) = 1;
    end
    if (Links.data(1,3)*d.z <= Init.z)
        alt_hist.lower(Links.data(1,3)+1) = 1;
    end
else
    alt_hist.data(Links.data(step,6)+1) = alt_hist.data(Links.data(step,6)+1) +1;
    if(Links.data(step,6)*d.z >= Init.z)
        alt_hist.upper(Links.data(step,6)+1) = alt_hist.upper(Links.data(step,6)+1) +1;
    end
    if (Links.data(step,6)*d.z <= Init.z)
        alt_hist.lower(Links.data(step,6)+1) = alt_hist.lower(Links.data(step,6)+1) +1;
    end
end
