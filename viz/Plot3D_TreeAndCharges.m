% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: Plot3D_TreeAndCharges.m                                     %
%    Purpose: Visualizes the discharge with the charge density overlaid.  %
%             Saves the initial and final structure. Can be turned into a %
%             video that tracks the propagation.                          %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: April 27, 2022                                              %
%    Updates:  November 2022 - Updated to account for complex charge      %
%                              density distributions.                     %
%              December 2024 - Integrated specifySimDetails.m function.   %
%             September 2025 - Renamed (formerly LightningVisual.m).      %
%               October 2025 - Integrated checkMagnitude.m function and   %
%                              tin-can BCs.                               %
%              December 2025 - Integrated the changes made to the         %
%                              specifySimDetails.m function and the newly %
%                              introduced setUpAxes.m function.           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close all
clearvars -except sims

%% Options to define variables that are typically input during run-time:
is.Rec            = 'N'; % Would you like to record the lightning propagation as a movie? (Y / N)
is.updateRho      = 'N'; % Would you like to update the charge distribution coloring for every saved step? (Y / N)
is.highResolution = 'N'; % Would you like to save the final image as a very high resolution image? (Y / N)
is.monoChrome     = 'N'; % Would you like the plot to be monochrome? (Y / N)

%% Fully-automated onwards:
if ~exist('sims','var')
    specifySimDetails;
end 

cd ../results

% Load data files
load('InitPoint.dat',        '-ascii');
rho.data  = load('rhoAmb.dat',           '-ascii');
Links.ID  = load('EstablishedLinks.dat', '-ascii');
stepsaves = abs(load('step3d.dat',       '-ascii'));
polarity  = load('TransportedRhoEnd.dat','-ascii');

if isempty(Links.ID)
    fprintf('\n*** Plot3D_TreeAndCharges.m cannot be executed with current EstablishedLinks.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot3D_TreeAndCharges.m script. ***\n');
end

% Derive main parameters (number of nodes, spacings, and domain size):
S.x = InitPoint(1);   % _m
S.y = InitPoint(2);   % _m
S.z = InitPoint(3);   % _m
S.R = InitPoint(4);   % _m

S.i = round(S.x/sims.domain.dx);
S.j = round(S.y/sims.domain.dy);
S.k = round(S.z/sims.domain.dz);

Links.Nb = size(Links.ID);
Links.Nb = Links.Nb(1);

clear dxyz
clear InitPoint
cd ../viz

rho.data = convertTo3d(rho.data,sims); % _C/_m^3

% Linear spaces for the three position dimensions:
x = (0:sims.domain.dx:sims.domain.maxx)'*sims.spatialFactor.Number;
y = (0:sims.domain.dy:sims.domain.maxy)'*sims.spatialFactor.Number;
z = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)'*sims.spatialFactor.Number;


[X,Y,Z]=meshgrid(x,y,z);
rho.max = max(max(max(rho.data)));
rho.min = min(min(min(rho.data)));

% Map ColorScale
color     = zeros(Links.Nb,3);
linestyle = repelem("-",Links.Nb);
if strcmp(is.monoChrome,'N')
    for ii=1:Links.Nb
        if polarity(ii)>0
            color(ii,:) = [1 0 0];
        elseif polarity(ii)<0
            color(ii,:) = [0 0 1];
        else
            color(ii,:) = [0 0 0];
            linestyle(ii) = "--";
        end
    end
end

% Set movie recording
if (strcmp(is.Rec,'Y') == 1)
    Movie = VideoWriter(strcat(sims.pathVideos,'/TreeAndCharges_',sims.objectName,'_',sims.objectType),'MPEG-4');
    if Links.Nb <= 60
        Movie.FrameRate = 2;
    elseif Links.Nb > 60 && Links.Nb <= 1000
        Movie.FrameRate = round(Links.Nb/30);
    else
        Movie.FrameRate = round(Links.Nb/60);
    end
    Movie.Quality = 100;
    open(Movie);
end

% Draw the tree
figure(1);
set(gcf,'Position',[0,0,1.2*sims.plotWidth,1.2*sims.plotHeight]);
set(gcf,'Resize','off')
hold on;

% Plots the cloud structure with the defined function below:
setUpAxes(sims,'xyz')
set(legend,'Position',[0.225 0.7 .5 .0375],'box','off')
set(gcf,'Resize','off')

% Initialize distance traveled for lightning links:
distance = 0;

% For-loop to plot lightning discharge links:
for ii=0:Links.Nb
    if ii == 0
        plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);
    else
        if strcmp(is.updateRho,'Y') == 1
            if mod((ii-1),stepsaves) == 0 || ii == Links.Nb
                allPatches = findall(gcf,'type','patch');
                delete(allPatches);
                if ii == Links.Nb
                    rho.data = load('../results/rho3d.dat','-ascii');
                else
                    rho.data = load(['../results/rho3d',num2str(ii-1),'.dat'],'-ascii');
                end
                rho.data = convertTo3d(rho.data,sims); % nC/_m^3
                if rho.max < max(max(max(rho.data))) || rho.max > 10*max(max(max(rho.data)))
                    rho.max = max(max(max(rho.data)));
                end
                if abs(rho.min) < abs(min(min(min(rho.data)))) || abs(rho.min) > 10*abs(min(min(min(rho.data))))
                    rho.min = min(min(min(rho.data)));
                end
                plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);
            end
        end    
        % Defining initial position of lightning link (m)
        x1 = Links.ID(ii,1)*sims.domain.dx;
        y1 = Links.ID(ii,2)*sims.domain.dy;
        z1 = Links.ID(ii,3)*sims.domain.dz+sims.domain.gnd; 
    
        % Defining final position of lightning link (m)
        x2 = Links.ID(ii,4)*sims.domain.dx;
        y2 = Links.ID(ii,5)*sims.domain.dy;
        z2 = Links.ID(ii,6)*sims.domain.dz+sims.domain.gnd;
    
        % Keeps track of initiation height and min/max propagation height:
        if ii == 1
            initHeight = z1;
            if z1 > z2
                maxHeight = z1;
                minHeight = z2;
            else
                maxHeight = z2;
                minHeight = z1;
            end
        else
            if z2 > maxHeight
                maxHeight = z2;
            end
            if z2 < minHeight
                minHeight = z2;
            end
        end
        
        % Summing lightning link to overall distance of link traveled (m):
        distance = distance + sqrt(((x2-x1)^2) + ((y2-y1)^2) + ((z2-z1)^2));
        
        % Plotting initiation point:
        if ii == 1
            scatter3(x1*sims.spatialFactor.Number,y1*sims.spatialFactor.Number,z1*sims.spatialFactor.Number,100,'filled','MarkerFaceColor','k','HandleVisibility','off','LineWidth',2);
        end

        % Plotting link:
        plot3([x1, x2]*sims.spatialFactor.Number,[y1, y2]*sims.spatialFactor.Number,[z1, z2]*sims.spatialFactor.Number,'Color',color(ii,:),'LineStyle',linestyle(ii),'HandleVisibility','off','LineWidth',2);
    end
    box on
    title(strcat(sims.objectType," discharge after ", int2str(ii) ," step(s)"),'FontSize',28,'FontWeight','bold','Interpreter','latex');
    if(strcmp(is.Rec,'Y') == 1)
        frame = getframe(gcf);
        writeVideo(Movie,frame);
        if ii == Links.Nb
            close(Movie);
        end
    end
end

% Outputs general information about the discharge's propagation:
fprintf(strcat('\n\t',sims.objectType," reaches minimum of ",num2str(minHeight*sims.spatialFactor.Number,'%.2f')," ",sims.spatialFactor.Prefix,'meters (',num2str((minHeight-initHeight)*sims.spatialFactor.Number,'%+.2f')," ",sims.spatialFactor.Prefix,'meters)\n'));
fprintf(strcat('\n\t',sims.objectType," initiated at ",num2str(initHeight*sims.spatialFactor.Number,'%.2f')," ",sims.spatialFactor.Prefix,'meters\n'));
fprintf(strcat('\n\t',sims.objectType," reaches maximum of ",num2str(maxHeight*sims.spatialFactor.Number,'%.2f')," ",sims.spatialFactor.Prefix,'meters (',num2str((maxHeight-initHeight)*sims.spatialFactor.Number,'%+.2f')," ",sims.spatialFactor.Prefix,'meters)\n'));
fprintf(strcat('\n\t',sims.objectType," has propagated a total of ",num2str(distance*sims.spatialFactor.Number,'%.2f')," ",sims.spatialFactor.Prefix,'meters\n'));

hold off;
%title('(b)','Interpreter','latex','FontSize',32,'Units','normalized');
%titleInfo = get(gca,'title');
%set(titleInfo,'Position', [((titleInfo.Extent(3))-(((titleInfo.Parent.InnerPosition(1)/titleInfo.Parent.InnerPosition(3)))/titleInfo.Parent.Position(3))-(titleInfo.Parent.Position(1)/titleInfo.Parent.OuterPosition(3))) 1.004 0]);
title(strcat("Simulated ", sims.objectType," Discharge: ",sims.objectName),'FontSize',28,'FontWeight','bold','Interpreter','latex');    

% If the 'export_fig' function is assigned to the pathtool:
if exist('export_fig') == 2 && strcmp(is.highResolution,'Y') == 1
    pause
    currentFolder = pwd;
    cd(sims.pathPNGs);
    export_fig HighRes_Discharge.png -transparent -m8
    cd(currentFolder);
else
    exportgraphics(gcf,strcat(sims.pathPNGs,'/Lightning_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
end