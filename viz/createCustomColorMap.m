% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: createCustomColorMap.m                                     %
%    Purpose: Creates a custom colormap for scalar or field values, where %
%             negative values are represented by blue and positive values %
%             are represented by red.                                     %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: February 22, 2022                                           %
%    Updates: April 2022 - Updated to account for desired transparency.   %
%             February 2025 - Introduced colormap for field values.       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function customColorMap = createCustomColorMap(valuetype,alphavalue)
    customColorMap = zeros([161,3]);
    fullValue = length(customColorMap(:,1));
    halfValue = ((fullValue-1)/2);
    if strcmp(valuetype,'scalar') == 1 || strcmp(valuetype,'Scalar') == 1
        % Version 1: blue to white to red
        customColorMap((halfValue+1),:) = [0.75 0.75 0.75];
        increasing_value = linspace((1-alphavalue),1,halfValue);
        decreasing_value = linspace(1,(1-alphavalue),halfValue);
        steady_value = ones([halfValue,1]);
        customColorMap(1:halfValue,1) = increasing_value;
        customColorMap(1:halfValue,2) = increasing_value;
        customColorMap(1:halfValue,3) = steady_value;
        customColorMap((halfValue+2):end,1) = steady_value;
        customColorMap((halfValue+2):end,2) = decreasing_value;
        customColorMap((halfValue+2):end,3) = decreasing_value;
    elseif strcmp(valuetype,'field') == 1 || strcmp(valuetype,'Field') == 1
        % Version 2: blue to yellow to red
        quarterValue = floor(halfValue/2)+1;
        steady_value = ones([halfValue,1]);
        steady_quarter = ones([quarterValue,1]);
        half_decrease = linspace(1,0,halfValue+1);
        quarter_decrease = linspace(1,0,quarterValue);
        quarter_increase = linspace(0,1,quarterValue);
        customColorMap(1:quarterValue,3) = quarter_decrease;
        customColorMap(1:quarterValue,2) = quarter_increase;
        customColorMap(quarterValue:halfValue+1,2) = steady_quarter;
        customColorMap(quarterValue:halfValue+1,1) = quarter_increase;
        customColorMap(halfValue+2:end,1) = steady_value;
        customColorMap(halfValue+1:end,2) = half_decrease;
    end
end