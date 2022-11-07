% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: createRedBlueColorMap.m                                      %
% Purpose: Creates a custom color map for negative charge densities (blue)%
%          and positive charge densities (red). Assumed preference for    %
%          colorblind-friendly mode, 'white'.                             %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: February 22, 2022                                           %
% Last Update: April 4, 2022 - Updated to account for an alpha value.     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function customColorMap = createRedBlueColorMap(neutralcolor,alphavalue)
    customColorMap = zeros([101,3]);
    %customColorMap = zeros([51,3]);
    customColorMap((((length(customColorMap(:,1))-1)/2)+1),:) = [0.75 0.75 0.75];
    if strcmp(neutralcolor,'white') == 1 || strcmp(neutralcolor,'White') == 1
        % Version 1: blue to white to red
        increasing_value = linspace((1-alphavalue),1,((length(customColorMap(:,1))-1)/2));
        decreasing_value = linspace(1,(1-alphavalue),((length(customColorMap(:,1))-1)/2));
        steady_value = ones([((length(customColorMap(:,1))-1)/2),1]);
        customColorMap(1:((length(customColorMap(:,1))-1)/2),1) = increasing_value;
        customColorMap(1:((length(customColorMap(:,1))-1)/2),2) = increasing_value;
        customColorMap(1:((length(customColorMap(:,1))-1)/2),3) = steady_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,1) = steady_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,2) = decreasing_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,3) = decreasing_value;
    else
        % Version 2: blue to grey to red
        increasing4zero_value = linspace(0,0.75,((length(customColorMap(:,1))-1)/2));
        decreasing4zero_value = linspace(0.75,0,((length(customColorMap(:,1))-1)/2));
        increasing4one_value = linspace(0.75,1,((length(customColorMap(:,1))-1)/2));
        decreasing4one_value = linspace(1,0.75,((length(customColorMap(:,1))-1)/2));
        customColorMap(1:((length(customColorMap(:,1))-1)/2),1) = increasing4zero_value;
        customColorMap(1:((length(customColorMap(:,1))-1)/2),2) = increasing4zero_value;
        customColorMap(1:((length(customColorMap(:,1))-1)/2),3) = decreasing4one_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,1) = increasing4one_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,2) = decreasing4zero_value;
        customColorMap((((length(customColorMap(:,1))-1)/2)+2):end,3) = decreasing4zero_value;
    end
end