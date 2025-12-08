% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: string2link.m                                               %
%    Purpose: Converts the string of N-characters that is returned by the %
%             characterizePropagation function into two N-by-3 arrays     %
%             to track the xyz start and finish points connecting the     %
%             initiation point to the final node. Essentially operates as %
%             the inverse of the characterizePropagation function for     %
%             plotting purposes. Originally embeded exclusively within    %
%             Plot1D_CurrentEstimate.m, but applications have since been  %
%             expanded.                                                   %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: November 29, 2025                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [starts, ends] = string2link(path, initiation)
    currentpath = convertStringsToChars(path);
    starts = zeros([length(currentpath),3]);
    ends = starts;
    starts(1,:) = initiation;
    for N = 1:1:length(currentpath)
        switch currentpath(N)
            case 'A'
                direction = [1 -1 1];
            case 'B'
                direction = [1 0 1];
            case 'C'
                direction = [1 1 1];
            case 'D'
                direction = [1 -1 0];
            case 'E'
                direction = [1 0 0];
            case 'F'
                direction = [1 1 0];
            case 'G'
                direction = [1 -1 -1];
            case 'H'
                direction = [1 0 -1];
            case 'I'
                direction = [1 1 -1];
            case 'J'
                direction = [0 -1 1];
            case 'K'
                direction = [0 0 1];
            case 'L'
                direction = [0 1 1];
            case 'M'
                direction = [0 -1 0];
            case 'N'
                direction = [0 1 0];
            case 'O'
                direction = [0 -1 -1];
            case 'P'
                direction = [0 0 -1];
            case 'Q'
                direction = [0 1 -1];
            case 'R'
                direction = [-1 -1 1];
            case 'S'
                direction = [-1 0 1];
            case 'T'
                direction = [-1 1 1];
            case 'U'
                direction = [-1 -1 0];
            case 'V'
                direction = [-1 0 0];
            case 'W'
                direction = [-1 1 0];
            case 'X'
                direction = [-1 -1 -1];
            case 'Y'
                direction = [-1 0 -1];
            case 'Z'
                direction = [-1 1 -1];
        end
        ends(N,:) = starts(N,:) + direction;
        if N~=length(currentpath)
            starts(N+1,:) = ends(N,:);
        end
    end 
end