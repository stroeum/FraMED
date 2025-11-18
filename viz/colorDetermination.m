% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: colorDetermination.m                                        %
%    Purpose: Determines the index and RGB triplet for a particular value %
%             on the maximum value for a quantity's range.                %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: February 22, 2022                                           %
%    Updates: October 2025 - Updated index range to match the values      %
%                            reflected in createCustomColorMap.m          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [colorIndex,colorRGB] = colorDetermination(currentValue,maxValue,colorMap)
    % Index of RGB triplet (range of 1-151, 51 is neutral/grounded):
    colorIndex = round((currentValue/abs(maxValue))*((length(colorMap(:,1))-1)/2))+(((length(colorMap(:,1))-1)/2)+1);
    % RGB triplet value:
    colorRGB = colorMap(colorIndex,:);
end