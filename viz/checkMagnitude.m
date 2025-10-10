% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: checkMagnitude.m                                            %
%    Purpose: Determines the appropriate unit naming convention based on  %
%             the magnitude of the numbers checked.                       %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: September 25, 2025                                          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function factor = checkMagnitude(number)
    % Determines the best order of magnitude based on an array of values:
    magnitude = log10(nonzeros(abs(number)));
    magBase3  = floor(magnitude/3);
    order = max(magBase3);
    
    % Limits possible values to the units available:
    if order > 3
        order = 3;
    elseif order < -5
        order = -5;
    end
    
    % Converts from SI to custom magnitude:
    factor.Number = 10^(-3*order);
    switch order
        case -5
            factor.Unit   = 'f';
            factor.Prefix = 'femto';
        case -4
            factor.Unit   = 'p';
            factor.Prefix = 'pico';
        case -3
            factor.Unit   = 'n';
            factor.Prefix = 'nano';
        case -2
            factor.Unit   = 'u';
            factor.Prefix = 'micro';
        case -1
            factor.Unit   = 'm';
            factor.Prefix = 'milli';
        case 0
            factor.Unit   = '';
            factor.Prefix = '';
        case 1
            factor.Unit   = 'k';
            factor.Prefix = 'kilo';
        case 2
            factor.Unit   = 'M';
            factor.Prefix = 'mega';
        case 3
            factor.Unit   = 'G';
            factor.Prefix = 'giga';
        otherwise
            return;
    end
    
    % Allows for appropriate rendering in the LaTeX environment:
    if order == -2
        factor.LaTeX = '$\mu$';
    else
        factor.LaTeX = factor.Unit;
    end
end