% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: plottingChargeRegions.m                                      %
% Purpose: Visualizes the charged cloud structure with custom color map.  %
%          Outputs a figure to the screen but does not save it to a file  %
%          automatically. LaTeX-ified April 29, 2022. Changed view angle  %
%          and several notations on February 20, 2024.                    %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: February 22, 2022                                           %
% Last Update: February 20, 2024                                          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function plottingChargeRegions(colorbarRange,alphaValue,rhoDataOG,Xval,Yval,Zval)
    %Creates a unique colormap to represent the charge regions:
    rgbValues = createRedBlueColorMap(colorbarRange,1);
    rgbValuesAdjusted = createRedBlueColorMap(colorbarRange,alphaValue);
   
    % Determines the range of the colorbar among other factors:
    tol = ceil(log10(round(max(max(max(abs(rhoDataOG.data)))),1,'significant')/(10^2)));    
    rhoData.data = round(rhoDataOG.data,-tol);
    uniqueRhos = unique(nonzeros(rhoData.data));
    rhoData.max = max(uniqueRhos);
    rhoData.min = min(uniqueRhos);

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
        [nullInd,~] = colorDetermination(trueUniqueRhos(i),max(abs(uniqueRhos)),rgbValues);
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
        [colorIndices(j),colorVertices(j,:)] = colorDetermination(uniqueRhos(j),max(abs(uniqueRhos)),rgbValues);
        
        % If the region is not neutrally charged:
        if colorIndices(j) ~= (((length(rgbValues(:,1))-1)/2)+1)
            % Is the current value a 'truly unique' value?
            [included, location] = ismember(uniqueRhos(j),trueUniqueRhos);
            % If so, add it to the legend:
            if included == 1
                %fprintf(['\n\n***included == 1***\nuniqueRhos(j) = ',num2str(uniqueRhos(j)),'\ntrueUniqueRhos(location) = ',num2str(trueUniqueRhos(location)),'\nuniqueRhos(location) = ',num2str(uniqueRhos(location))]);
                if trueUniqueRhos(location)>0
                    p1 = patch(isosurface(Xval,Yval,Zval,-1*pagetranspose(rhoData.data),-trueUniqueRhos(location)));
                else
                    p1 = patch(isosurface(Xval,Yval,Zval,pagetranspose(rhoData.data),trueUniqueRhos(location)));
                end
                if length(trueUniqueRhos)>3
                    if colorIndices(j) <= (middle+midRange) && colorIndices(j) >= (middle-midRange)
                        chosenAlpha = alphaValue/3;
                    else
                        chosenAlpha = alphaValue;
                    end
                    if uniqueRhos(location) == max_rho_value && min_rho_value > 0
                        set(p1,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','on','DisplayName',['$$0 \leq \rho \leq$$ ',num2str(trueUniqueRhos(location)),' nC/m$^3$']);
                        isonormals(Xval,Yval,Zval,pagetranspose(rhoData.data),p1);
                        drawnow
                    elseif uniqueRhos(location) == min_rho_value && max_rho_value < 0
                        set(p1,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','on','DisplayName',['$$0 \geq \rho \geq$$ ',num2str(trueUniqueRhos(location)),' nC/m$^3$']);
                        isonormals(Xval,Yval,Zval,pagetranspose(rhoData.data),p1);
                        drawnow
                    else
                        set(p1,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','off');
                        isonormals(Xval,Yval,Zval,pagetranspose(rhoData.data),p1);
                        drawnow
                    end
                else
                    set(p1,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',alphaValue,'HandleVisibility','on','DisplayName',['Charge Density $$\approx$$ ',num2str(trueUniqueRhos(location)),' nC/m$^3$']);
                    isonormals(Xval,Yval,Zval,pagetranspose(rhoData.data),p1);
                    drawnow
                end
            else
                %patch(isosurface(Xval,Yval,Zval,pagetranspose(rhoData.data),uniqueRhos(j)),'FaceColor',colorVertices(j,:),'EdgeColor',colorVertices(j,:),'FaceAlpha',alphaValue,'HandleVisibility','off'); % set the color, mesh and transparency level of the surface
            end
        end
    end
    
    % Resolves formatting issues with the custom colorbar:
    cmap = (rgbValues+rgbValuesAdjusted+rgbValuesAdjusted)/3;
    caxis([-max(abs(uniqueRhos)) max(abs(uniqueRhos))]);
    colormap(cmap);
    colorbar;
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
end