% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: processCurrents.m                                           %
%    Purpose: Calculates the current through every link of the structure  %
%             at all steps.                                               %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [links,nodes] = processCurrents(links,nodes)
    for N = 1:1:links.total
        % Processes and calculates the current through each link at the respective location:
        endingNode = N + 1;
        startingNode = links.connections.nodes(N) + 1;
        while links.connections.num(N,startingNode) + links.connections.num(N,endingNode) > 2
            [values,indices]=sort(links.connections.num(N,1:N));
            focus = indices(values==1);
            for ii = 1:1:length(focus)
                if (links.connections.num(N,startingNode) > 1)
                    singleNode = focus(ii);
                    if singleNode > 1 && links.connections.num(N,links.connections.nodes(singleNode-1)+1)>=1
                        connectedNode = links.connections.nodes(singleNode-1)+1;
                    else
                        [tempvals,tempind] = sort(links.connections.nodes(1:N)');
                        subfocus = tempind(tempvals==(singleNode-1));
                        for jj = 1:1:length(subfocus)
                            if links.connections.num(N,subfocus(jj)+1)>=1
                                connectedNode = subfocus(jj)+1;
                            end
                        end
                    end
                    links.connections.deltas(N,singleNode) = links.connections.deltas(N,singleNode) + nodes.rho.deltas(N,singleNode);
                    links.connections.deltas(N,connectedNode) = links.connections.deltas(N,connectedNode) + links.connections.deltas(N,singleNode);
                    links.connections.num(N,singleNode) = links.connections.num(N,singleNode) - 1;
                    links.connections.num(N,connectedNode) = links.connections.num(N,connectedNode) - 1;
    
                    if links.connections.num(N,singleNode) == 0
                        links.connections.currents(N,(max([singleNode; connectedNode])-1)) = links.connections.deltas(N,singleNode);
                    end
                end
            end
        end
        links.connections.deltas(N,startingNode) = links.connections.deltas(N,startingNode) + nodes.rho.deltas(N,startingNode);
        links.connections.deltas(N,endingNode) = links.connections.deltas(N,startingNode) + nodes.rho.deltas(N,endingNode);
        links.connections.num(N,startingNode) = links.connections.num(N,startingNode) - 1;
        links.connections.num(N,endingNode) = links.connections.num(N,endingNode) - 1;
        links.connections.currents(N,N) = links.connections.deltas(N,startingNode);
    end
end