function plottingChargeRegions(colorbarRange,alphaValue,rhoData,Xval,Yval,Zval,xval,yval,zval)
    %Creates a unique colormap to represent the charge regions:
    rgbValues = createRedBlueColorMap(colorbarRange); 
    % Determines the range of the colorbar among other factors:
    uniqueRhos = nonzeros(unique(rhoData.data));
    colorIndices = zeros([length(uniqueRhos),1]);
    colorVertices = zeros([length(uniqueRhos),3]);
    % Determines truly unique values for legend usage:
    trueUniqueRhos = nonzeros(uniquetol(rhoData.data,10^-4));
    % Removes values that are approximately zero (i.e. neutral):
    [null,nullLoc] = ismember(0,round(trueUniqueRhos,4));
    if null == 1
        trueUniqueRhos(nullLoc)=[];
    end
    Legend = cell(length(trueUniqueRhos),1);

    % Plots isocharge regions in a representative color:
    for j = 1:1:length(uniqueRhos)
        [colorIndices(j),colorVertices(j,:)] = colorDetermination(uniqueRhos(j),max(abs(uniqueRhos)),rgbValues);
        % If the region is not neutrally charged:
        if colorIndices(j) ~= 51
            % Is the current value a 'truly unique' value?
            [included, location] = ismember(uniqueRhos(j),trueUniqueRhos);
            % If so, add it to the legend:
            if included == 1
                uniquePatch = patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(j)),'FaceAlpha',alphaValue,'FaceColor',colorVertices(j,:),'EdgeColor','none'); 
                Legend{location}=['Charge Density \approx ',num2str(trueUniqueRhos(location)),' nC/m^3'];
                isonormals(xval,yval,zval,rhoData.data,uniquePatch);
                set(uniquePatch,'FaceAlpha',0.06)% set the color, mesh and transparency level of the surface
            else
                isonormals(xval,yval,zval,rhoData.data,patch(isosurface(Xval,Yval,Zval,rhoData.data,uniqueRhos(j)),'FaceColor',colorVertices(j,:),'EdgeColor','none','FaceAlpha',alphaValue,'HandleVisibility','off')); % set the color, mesh and transparency level of the surface
            end
        end
    end
    % Resolves formatting issues with the colorbar:
    cmap = rgbValues;
    caxis([-max(abs(uniqueRhos)) max(abs(uniqueRhos))]);
    colormap(cmap);
    colorbar;
    % Adds the legend to the plot:
    [~,h_legend] = legend(Legend,'Box','off');
    PatchInLegend = findobj(h_legend, 'type', 'patch');
    set(PatchInLegend(:), 'FaceAlpha', 1);
    % Additional formatting aspects:
    view([33 10])
    ax = gca;
    ax.Colorbar.Title.String = "nC/m^3";
    ax.Color;
    xlabel('x-position (km)');
    ylabel('y-position (km)');
    zlabel('z-position (km)');
    grid on
end