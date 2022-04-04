% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: plottingChargeRegions.m                                      %
% Purpose: Visualizes the charged cloud structure with custom color map.  %
%          Outputs a figure to the screen but does not save it to a file  %
%          automatically.                                                 %
% Author: Annelisa Esparza                                                %
% Contact: aesparza2014@my.fit.edu                                        %
% Added Date: February 22, 2022                                           %
% Last Update: April 4, 2022                                              %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function plottingChargeRegions(colorbarRange,alphaValue,rhoDataOG,Xval,Yval,Zval,xval,yval,zval)
    %Creates a unique colormap to represent the charge regions:
    rgbValues = createRedBlueColorMap(colorbarRange); 
   
    % Determines the range of the colorbar among other factors:
    tol = ceil(log10(round(max(max(max(abs(rhoDataOG.data)))),1,'significant')/10^4));    
    rhoData.data = round(rhoDataOG.data,-tol);
    uniqueRhos = unique(nonzeros(rhoData.data));
    rhoData.max = max(uniqueRhos);
    rhoData.min = min(uniqueRhos);

    % Prevents misread of charge density:
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
    for testloop2 = 1:1:length(uniqueRhos)
        if testingvalues(testloop2)==1
            rhoData.data(:,:,testinglocation(testloop2)-1)=rhoData.data(:,:,testinglocation(testloop2));
        end
    end

    colorIndices = zeros([length(uniqueRhos),1]);
    colorVertices = zeros([length(uniqueRhos),3]);
    [colorIndices(1),colorVertices(1,:)] = colorDetermination(uniqueRhos(1),max(abs(uniqueRhos)),rgbValues);
    [colorIndices(2),colorVertices(2,:)] = colorDetermination(uniqueRhos(2),max(abs(uniqueRhos)),rgbValues);
    % Determines truly unique values for legend usage:
    trueUniqueRhos = uniqueRhos;
    % Removes values that are approximately zero (i.e. neutral):
    for i = length(trueUniqueRhos):-1:1
        [nullInd,~] = colorDetermination(trueUniqueRhos(i),max(abs(uniqueRhos)),rgbValues);
        if nullInd == 51
            trueUniqueRhos(i)=[];
        end
    end
    view([33 10])
    % Plots isocharge regions in a representative color:
    for j = length(uniqueRhos):-1:1
        [colorIndices(j),colorVertices(j,:)] = colorDetermination(uniqueRhos(j),max(abs(uniqueRhos)),rgbValues);
        % If the region is not neutrally charged:
        if colorIndices(j) ~= 51
            % Is the current value a 'truly unique' value?
            [included, location] = ismember(uniqueRhos(j),trueUniqueRhos);
            % If so, add it to the legend:
            if included == 1
                if location == 1
                    disp(rhoDataOG.min>uniqueRhos(location))
                    if rhoDataOG.max<uniqueRhos(location)
                        p1 = patch(isosurface(Xval,Yval,Zval,-1*rhoData.data,-uniqueRhos(location)));
                    else
                        p1 = patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(location)));
                    end
                    set(p1,'FaceColor',colorVertices(location,:),'EdgeAlpha',0,'FaceAlpha',0.75,'HandleVisibility','on','DisplayName',['Charge Density \approx ',num2str(trueUniqueRhos(location)),' nC/m^3']);
                    isonormals(Xval,Yval,Zval,rhoData.data,p1);
                    drawnow
                elseif location == 2
                    if rhoDataOG.max<uniqueRhos(location)
                        p2 = patch(isosurface(Xval,Yval,Zval,-1*rhoData.data,-uniqueRhos(location)));
                    elseif rhoDataOG.min>uniqueRhos(location)
                        p2 = patch(isosurface(Xval,Yval,Zval,-1*rhoData.data,-uniqueRhos(location)));
                    else
                        p2 = patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(location)));
                    end
                    set(p2,'FaceColor',colorVertices(location,:),'EdgeAlpha',0,'FaceAlpha',0.75,'HandleVisibility','on','DisplayName',['Charge Density \approx ',num2str(trueUniqueRhos(location)),' nC/m^3']);
                    isonormals(Xval,Yval,Zval,rhoData.data,p2);
                    drawnow
                elseif location == 3
                    p3 = patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(j)));
                    set(p3,'FaceColor',colorVertices(j,:),'EdgeColor',colorVertices(j,:),'FaceAlpha',0.75,'HandleVisibility','on','DisplayName',['Charge Density \approx ',num2str(trueUniqueRhos(location)),' nC/m^3']);
                    isonormals(Xval,Yval,Zval,rhoData.data,p3);
                    drawnow
                end
            else
                patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(j)),'FaceColor',colorVertices(j,:),'EdgeColor',colorVertices(j,:),'FaceAlpha',0.5,'HandleVisibility','off'); % set the color, mesh and transparency level of the surface
            end
        end
    end
    
    % Resolves formatting issues with the colorbar:
    cmap = rgbValues;
    caxis([-max(abs(uniqueRhos)) max(abs(uniqueRhos))]);
    colormap(cmap);
    colorbar;
    % Adds the legend to the plot:
    legend('boxoff');
    ax = gca;
    ax.Colorbar.Title.String = "nC/m^3";
    ax.Color;
    xlabel('x-position (km)');
    ylabel('y-position (km)');
    zlabel('z-position (km)');
    grid on
end