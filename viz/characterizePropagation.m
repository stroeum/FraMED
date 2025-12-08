% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: characterizePropagation.m                                   %
%    Purpose: Identifies the link's position with respect to the          %
%             initiation point. Originally embeded exclusively within     %
%             Plot1D_CurrentEstimate.m, but applications have since been  %
%             expanded.                                                   %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: November 29, 2025                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function letter = characterizePropagation(links, N)
    % Letters A - I:
    if links.end(N,1) == (links.start(N,1)+1)
        % Letters A - C:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'A';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'B';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'C';
            end
        % Letters D - F:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'D';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'E';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'F';
            end
        % Letters G - I:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'G';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'H';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'I';
            end
        end
    % Letters J - Q:
    elseif links.end(N,1) == links.start(N,1)
        % Letters J - L:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'J';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'K';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'L';
            end
        % Letters M - N:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'M';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'N';
            end
        % Letters O - Q:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'O';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'P';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'Q';
            end
        end
    % Letters R - Z:
    elseif links.end(N,1) == (links.start(N,1)-1)
        % Letters R - T:
        if links.end(N,3) == (links.start(N,3)+1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'R';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'S';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'T';
            end
        % Letters U - W:
        elseif links.end(N,3) == links.start(N,3)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'U';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'V';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'W';
            end
        % Letters X - Z:
        elseif links.end(N,3) == (links.start(N,3)-1)
            if links.end(N,2) == (links.start(N,2)-1)
                letter = 'X';
            elseif links.end(N,2) == links.start(N,2)
                letter = 'Y';
            elseif links.end(N,2) == (links.start(N,2)+1)
                letter = 'Z';
            end
        end
    end
end