% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: colorDetermination.m                                         %
% Purpose: Determines the index and RGB triplet for a particular value    %
%          on the maximum value for a quantity's range.                   %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: February 22, 2022                                           %
% Last Update: N/A                                                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [colorIndex,colorRGB] = colorDetermination(currentValue,maxValue,colorMap)
    % Index of RGB triplet (range of 1-101, 51 is neutral/grounded):
    colorIndex = round((currentValue/abs(maxValue))*50)+51;
    % RGB triplet value:
    colorRGB = colorMap(colorIndex,:);
end