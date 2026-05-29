% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: processCurrents.m                                           %
%    Purpose: Calculates the current through every link of the structure  %
%             at all steps.                                               %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [links,rhos] = processCurrents(links,rhos)
    links.direction = NaN([links.total,links.total]);
    for ll = 1:1:links.total
        % Processes and calculates the current through each link at the respective location:
        endingNode = ll + 1;
        startingNode = links.connections.nodes(ll) + 1;
        % While there are nodes 
        while links.connections.num(ll,startingNode) + links.connections.num(ll,endingNode) > 2
            % Sort the number of connections from low to high:
            [values,indices]=sort(links.connections.num(ll,1:ll));
            % Isolate those that have only one connection:
            focus = indices(values==1);
            for ii = 1:1:length(focus)
                if (links.connections.num(ll,startingNode) > 1)
                    singleNode = focus(ii);
                    if singleNode > 1 && links.connections.num(ll,links.connections.nodes(singleNode-1)+1)>1
                        connectedNode = links.connections.nodes(singleNode-1)+1;
                    else
                        [tempvals,tempind] = sort(links.connections.nodes(1:ll)');
                        subfocus = tempind(tempvals==(singleNode-1));
                        for jj = 1:1:length(subfocus)
                            if links.connections.num(ll,subfocus(jj)+1)>=1
                                connectedNode = subfocus(jj)+1;
                            end
                        end
                    end
                    links.connections.deltas(ll,singleNode) = links.connections.deltas(ll,singleNode) + rhos.deltas(ll,singleNode);
                    links.connections.deltas(ll,connectedNode) = links.connections.deltas(ll,connectedNode) + links.connections.deltas(ll,singleNode);
                    links.connections.num(ll,singleNode) = links.connections.num(ll,singleNode) - 1;
                    links.connections.num(ll,connectedNode) = links.connections.num(ll,connectedNode) - 1;
                    
                    if links.connections.num(ll,singleNode) == 0
                        if singleNode>1 && connectedNode == (links.connections.nodes(singleNode-1)+1)
                            links.direction(ll,singleNode-1) = -1;
                        elseif singleNode>1 && singleNode == (links.connections.nodes(connectedNode-1)+1)
                            links.direction(ll,connectedNode-1) = 1;
                        elseif singleNode == 1 && singleNode ~= (links.connections.nodes(connectedNode-1)+1)
                            links.direction(ll,connectedNode-1) = -1;
                        elseif singleNode == 1 && singleNode == (links.connections.nodes(connectedNode-1)+1)
                            links.direction(ll,connectedNode-1) = 1;
                        else
                            links.direction(ll,singleNode-1) = 1;
                        end
                        links.connections.currents(ll,(max([singleNode; connectedNode])-1)) = links.connections.deltas(ll,singleNode);
                    end
                end
            end
        end
        if startingNode>1 && singleNode == (links.connections.nodes(startingNode-1)+1)
            links.direction(ll,startingNode-1) = 1;
        end
        links.direction(ll,endingNode-1) = 1;
        links.connections.deltas(ll,startingNode) = links.connections.deltas(ll,startingNode) + rhos.deltas(ll,startingNode);
        links.connections.deltas(ll,endingNode) = links.connections.deltas(ll,startingNode) + rhos.deltas(ll,endingNode);
        links.connections.num(ll,startingNode) = links.connections.num(ll,startingNode) - 1;
        links.connections.num(ll,endingNode) = links.connections.num(ll,endingNode) - 1;
        links.connections.currents(ll,ll) = links.connections.deltas(ll,startingNode);
    end
end