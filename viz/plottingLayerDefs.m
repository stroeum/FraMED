% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   File Name: plottingLayerDefs.m                                        %
%     Purpose: Visualizes the charge layer with FraMED definitions.       %
%              Outputs a figure to the screen but does not save it to a   %
%              file automatically.                                        %
%      Author: Annelisa Esparza                                           %
%     Contact: annelisa.esparza@my.erau.edu                               %
%  Added Date: April 29, 2022                                             %
% Last Update: February 2025 - Updated to match the options available for %
%                              the createCustomColorMap function.        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function plottingLayerDefs(colorbarRange,alphaValue,rhoDataOG,Xval,Yval,Zval,R1,h1,Q1,pos1,R2,h2,Q2,pos2)
    % Ensures charge regions are properly labeled
    if Q1 < Q2
        tempQ = Q1;
        tempR = R1;
        tempH = h1;
        tempPos = pos1;
        Q1 = Q2;
        R1 = R2;
        h1 = h2;
        pos1 = pos2;
        Q2 = tempQ;
        R2 = tempR;
        h2 = tempH;
        pos2 = tempPos;
        clear tempQ tempR tempH tempPos
    end

    %Creates a unique colormap to represent the charge regions:
    rgbValues = createCustomColorMap(colorbarRange,1);
    rgbValuesAdjusted = createCustomColorMap(colorbarRange,alphaValue);
   
    % Determines the range of the colorbar among other factors:
    tol = ceil(log10(round(max(max(max(abs(rhoDataOG.data)))),1,'significant')/(10^2)));    
    rhoData.data = round(rhoDataOG.data,-tol);
    uniqueRhos = unique(nonzeros(rhoData.data));
    rhoData.max = max([rhoDataOG.max max(uniqueRhos)]);
    rhoData.min = min([rhoDataOG.min min(uniqueRhos)]);

    % Prevents misread of charge density when cloud height is less than dz:
    testingvalues = zeros([length(uniqueRhos),1]);
    testinglocation = zeros([length(uniqueRhos),1]);
    for testloop = 1:1:length(rhoData.data(1,1,:))
        testingunique = unique(nonzeros(rhoData.data(:,:,testloop)));
        if testingunique~=0
            [~, testloc] = ismember(testingunique,uniqueRhos);
            testingvalues(testloc) = testingvalues(testloc)+1;
            testinglocation(testloc) = testloop;
        end
    end

    % Copies data to a nearby height without overwriting pre-existing data:
    for testloop2 = 1:1:length(uniqueRhos)
        if testingvalues(testloop2)==1
            if max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)-1))))==0 && max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)+1))))==0
                rhoData.data(:,:,testinglocation(testloop2)-1)=rhoData.data(:,:,testinglocation(testloop2));
            elseif max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)-1))))~=0 && max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)+1))))==0
                rhoData.data(:,:,testinglocation(testloop2)+1)=rhoData.data(:,:,testinglocation(testloop2));
            elseif max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)-1))))==0 && max(abs(unique(rhoData.data(:,:,testinglocation(testloop2)+1))))~=0
                rhoData.data(:,:,testinglocation(testloop2)-1)=rhoData.data(:,:,testinglocation(testloop2));
            end
        end
    end
    % Determines truly unique values for legend usage:
    trueUniqueRhos = uniqueRhos;
    middle = (((length(rgbValues(:,1))-1)/2)+1);
    midRange = floor(((length(rgbValues(:,1))-1)/10))-1;
    % Removes values that are approximately zero (i.e. neutral):
    for i = length(trueUniqueRhos):-1:1
        [nullInd,~] = colorDetermination(trueUniqueRhos(i),max(abs([rhoData.min rhoData.max])),rgbValues);
        %if nullInd == (((length(rgbValues(:,1))-1)/2)+1)
        if nullInd <= (middle+round(midRange/2)) && nullInd >= (middle-round(midRange/2))
            trueUniqueRhos(i)=[];
        end
    end
    max_rho_value = max(trueUniqueRhos);
    min_rho_value = min(trueUniqueRhos);
    %fprintf(['Maximum density is ',num2str(max(trueUniqueRhos)),'\nMinimum density is ',num2str(min(trueUniqueRhos)),'\nLength of unique densities is ',num2str(length(trueUniqueRhos))]);
    % Sets 
    view(-45,9);
    % Plots isocharge regions in a representative color:
    colorIndices = zeros([length(uniqueRhos),1]);
    colorVertices = zeros([length(uniqueRhos),3]);
    for j = length(uniqueRhos):-1:1
        [colorIndices(j),colorVertices(j,:)] = colorDetermination(uniqueRhos(j),max(abs([rhoData.min rhoData.max])),rgbValues);
        
        % If the region is not neutrally charged:
        if colorIndices(j) ~= (((length(rgbValues(:,1))-1)/2)+1)
            % Is the current value a 'truly unique' value?
            [included, location] = ismember(uniqueRhos(j),trueUniqueRhos);
            % If so, add it to the legend:
            % If so, add it to the legend:
            if included == 1
                if uniqueRhos(location)>0
                    p1 = patch(isosurface(Xval,Yval,Zval,-1*rhoData.data,-uniqueRhos(location)));
                    set(p1,'FaceColor',colorVertices(location,:),'EdgeAlpha',0,'FaceAlpha',alphaValue,'HandleVisibility','on','DisplayName',[newline 'Radius: $R_1$ = ',num2str(R1),' km' newline 'Height: $h_1$ = ',num2str(h1),' km' newline 'Charge: $Q_1$ = ',num2str(Q1),' C' newline 'Center: $\left(x_1,y_1,z_1\right)$ = (',num2str(pos1(1)),' km, ',num2str(pos1(2)),' km, ',num2str(pos1(3)),' km)' newline]);
                else
                    p1 = patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(location)));
                    set(p1,'FaceColor',colorVertices(location,:),'EdgeAlpha',0,'FaceAlpha',alphaValue,'HandleVisibility','on','DisplayName',[newline 'Radius: $R_2$ = ',num2str(R2),' km' newline 'Height: $h_2$ = ',num2str(h2),' km' newline 'Charge: $Q_2$ = ',num2str(Q2),' C' newline 'Center: $\left(x_2,y_2,z_2\right)$ = (',num2str(pos2(1)),' km, ',num2str(pos2(2)),' km, ',num2str(pos2(3)),' km)' newline]);
                end
                isonormals(Xval,Yval,Zval,rhoData.data,p1);
                drawnow
            else
                patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(j)),'FaceColor',colorVertices(j,:),'EdgeColor',colorVertices(j,:),'FaceAlpha',alphaValue,'HandleVisibility','off'); % set the color, mesh and transparency level of the surface
            end
        end
    end
    
    % Resolves formatting issues with the custom colorbar:
    cmap = (rgbValues+rgbValuesAdjusted+rgbValuesAdjusted)/3;
    clim([-max(abs([rhoData.min rhoData.max max(abs(uniqueRhos))])) max(abs([rhoData.min rhoData.max max(abs(uniqueRhos))]))]);
    colormap(cmap);
    c = colorbar;
    c.Label.String = 'Charge Density (nC/m$^3$)';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 24;
    c.TickLabelInterpreter = 'latex';
    
    % Adds the legend to the plot:
    [~,h_legend]=legend('Location','north','Box','off','FontSize',18,'Interpreter','latex');
    PatchInLegend = findobj(h_legend,'type','patch');
    set(PatchInLegend(:),'FaceAlpha',((1-alphaValue)*alphaValue)+alphaValue);
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 20;
    xlabel('$x$-position (km)','Interpreter','latex','FontSize',22,'HorizontalAlignment','left');
    ylabel('$y$-position (km)','Interpreter','latex','FontSize',22,'HorizontalAlignment','right');
    zlabel('$z$-position (km)','Interpreter','latex','FontSize',24);
    grid on
    view(-45,9)
    if length(trueUniqueRhos)<=3
        pause
    end
end