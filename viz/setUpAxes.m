% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: setUpAxes.m                                                 %
%    Purpose: Consolidates the formatting for figures. Currently only     %
%             applicable for 3D visualizations.                           %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function setUpAxes(sims,variant,backgroundColor)
    if (~exist('backgroundColor', 'var'))
        backgroundColor = 'w';
    end
    switch variant
        case 'xyz'
            if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                gnd.color = [.75 .75 .75];
                P.x = [sims.domain.maxx sims.domain.minx sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number;
                P.y = [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number;
                P.z = [sims.domain.gnd sims.domain.gnd sims.domain.gnd sims.domain.gnd]*sims.spatialFactor.Number;
                patch(P.x, P.y, P.z, sims.domain.gnd*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off');
                if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                    patch(P.x, P.y, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    if strcmp(sims.BCtype,'TIN_CAN')
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.maxx sims.domain.maxx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.minx sims.domain.minx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.maxy sims.domain.maxy]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    end
                end
            end
            axis equal
            set(gca,'Color',backgroundColor,'TickLabelInterpreter', 'latex','FontSize',16)
            xlabel(strcat('$x$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',22,'HorizontalAlignment','left');
            ylabel(strcat('$y$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',22,'HorizontalAlignment','right');
            zlabel(strcat('$z$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',24);
            xlim([sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number);
            ylim([sims.domain.miny sims.domain.maxy]*sims.spatialFactor.Number);
            zlim([sims.domain.minz sims.domain.maxz]*sims.spatialFactor.Number);
            grid on
            view(-45,5)
        case 'xyz_hires'
            if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                gnd.color = [.75 .75 .75];
                P.x = [sims.domain.maxx sims.domain.minx sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number;
                P.y = [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number;
                P.z = [sims.domain.gnd sims.domain.gnd sims.domain.gnd sims.domain.gnd]*sims.spatialFactor.Number;
                patch(P.x, P.y, P.z, sims.domain.gnd*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off');
                if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                    patch(P.x, P.y, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    if strcmp(sims.BCtype,'TIN_CAN')
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.maxx sims.domain.maxx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.minx sims.domain.minx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.maxy sims.domain.maxy]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    end
                end
            end
            axis equal
            set(gca,'Color',backgroundColor,'TickLabelInterpreter', 'latex','FontSize',2*16)
            xlabel(strcat('$x$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',2*22,'HorizontalAlignment','left');
            ylabel(strcat('$y$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',2*22,'HorizontalAlignment','right');
            zlabel(strcat('$z$-position (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',2*24);
            xlim([sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number);
            ylim([sims.domain.miny sims.domain.maxy]*sims.spatialFactor.Number);
            zlim([sims.domain.minz sims.domain.maxz]*sims.spatialFactor.Number);
            grid on
            view(-45,5)
        case '3d'
            if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                gnd.color = [.75 .75 .75];
                P.x = [sims.domain.maxx sims.domain.minx sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number;
                P.y = [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number;
                P.z = [sims.domain.gnd sims.domain.gnd sims.domain.gnd sims.domain.gnd]*sims.spatialFactor.Number;
                patch(P.x, P.y, P.z, sims.domain.gnd*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off');
                if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                    patch(P.x, P.y, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number, [sims.domain.maxz sims.domain.maxz sims.domain.maxz sims.domain.maxz]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    if strcmp(sims.BCtype,'TIN_CAN')
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.maxx sims.domain.maxx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.miny sims.domain.miny]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                        patch([sims.domain.maxx sims.domain.maxx sims.domain.minx sims.domain.minx]*sims.spatialFactor.Number, [sims.domain.maxy sims.domain.maxy sims.domain.maxy sims.domain.maxy]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number, [sims.domain.gnd sims.domain.maxz sims.domain.maxz sims.domain.gnd]*sims.spatialFactor.Number,'FaceColor',gnd.color,'HandleVisibility','off'); 
                    end
                end
            end
            axis equal
            xlim([sims.domain.minx sims.domain.maxx]*sims.spatialFactor.Number);
            ylim([sims.domain.miny sims.domain.maxy]*sims.spatialFactor.Number);
            zlim([sims.domain.minz sims.domain.maxz]*sims.spatialFactor.Number);
            grid on
            view(-45,5)
    end
end