% Determines the index and RGB triplet for a particular value based on the
% maximum value for a quantity's range:
function [colorIndex,colorRGB] = colorDetermination(currentValue,maxValue,colorMap)
    % Index of RGB triplet (range of 1-101, 51 is neutral/grounded):
    colorIndex = round((currentValue/abs(maxValue))*50)+51;
    % RGB triplet value:
    colorRGB = colorMap(colorIndex,:);
end