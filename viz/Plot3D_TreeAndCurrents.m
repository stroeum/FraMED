% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: Plot3D_TreeAndCurrents.m                                    %
%    Purpose: Approximates the current through all links at all steps in  %
%             the structure with the assumption of a known, constant      % 
%             propagation speed defined by whether the discharge is a     %
%             streamer or leader. Outputs to the screen relevant values   %
%             and statistics similar to Plot1D_CurrentEstimate.m.         %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 8, 2025                                            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close all
clearvars -except sims

%% Options to define variables that are typically input during run-time:
is.Rec            = 'Y'; % Would you like to record the lightning propagation as a movie? (Y / N)
is.highResolution = 'N'; % Would you like to save the final image as a very high resolution image? (Y / N)

is.currentVisual   = 'log'; % Sets the scaling approach for the link's current intensity. Options: 'log' or 'linear'.
is.backgroundColor = 'w';   % Sets the plot's background color.
is.textColor       = 'k';   % Sets the plot's text color.

%% Fully-automated onwards:
fprintf('\n*** Executing Plot3D_TreeAndCurrents.m script. ***\n');  
if ~exist('sims','var')
    specifySimDetails;
end

% Read-in data for simulation and process the current at the tip:
[current, polarity, tau, links, nodes] = processTipCurrent(sims);

% Calculate the current through all links at all steps: 
[links, nodes] = processCurrents(links,nodes);

% Determine the scale:
temporalFactor    = checkMagnitude([links.timescale.min; links.timescale.max]);
currentFactor     = checkMagnitude(current.case(:));
minCurrentFactor  = checkMagnitude(current.min.val);
conversionFactor  = (10^(-9))*sims.domain.dx*sims.domain.dy*sims.domain.dz;
chargeTransported = conversionFactor.*abs(current.polarity(1:links.total));

% Output to the screen the results:
fprintf(strcat("\n\t",sims.objectType," (positive) propagation speed approximated as ",num2str(sims.vprop.pos/1000,'%.0f')," km/s."))
fprintf(strcat("\n\t",sims.objectType," (negative) propagation speed approximated as ",num2str(sims.vprop.neg/1000,'%.0f')," km/s."))
fprintf(strcat("\n\tDischarge timescale estimated to be between ",num2str(temporalFactor.Number*links.timescale.min,'%.1f')," - ",num2str(temporalFactor.Number*links.timescale.max,'%.1f')," ",temporalFactor.Unit,"s."))
fprintf(strcat("\n\tLongest branch traveled ",num2str(links.maxdistance*sims.spatialFactor.Number,'%.1f')," ",sims.spatialFactor.Unit,"m from the initiation point."))
fprintf(strcat("\n\tAverage currents (based on charge transfer):"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.avg.val,5),"\t",currentFactor.Unit,"A\t(mean, including all links,  sigma = ",num2str(currentFactor.Number*current.avg.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.partial.val,5),"\t",currentFactor.Unit,"A\t(mean, excluding final link, sigma = ",num2str(currentFactor.Number*current.partial.std,5)," ",currentFactor.Unit,"A)"));
fprintf(strcat("\n\t\t",num2str(currentFactor.Number*current.max.val,5),"\t",currentFactor.Unit,"A\t(maximum, link ",num2str(links.steps(current.max.index)),")"));
fprintf(strcat("\n\t\t",num2str(minCurrentFactor.Number*current.min.val,5),"\t",minCurrentFactor.Unit,"A\t(minimum, link ",num2str(links.steps(current.min.index)),")\n"));
if links.end(end,3)==0
    fprintf(strcat("\t\t",num2str(0.001*chargeTransported(end)./(links.pathdistance(end)/(polarity.neg(end)*sims.vprop.neg + polarity.pos(end)*sims.vprop.pos)),5),"\tkA\t(initiation-to-ground path)\n"))
end

%% Plotting the current density:
% Set movie recording
if (strcmp(is.Rec,'Y') == 1)
    Movie = VideoWriter(strcat(sims.pathVideos,'/',sims.objectName,'_',sims.objectType,'_CurrentMovement_Video'),'MPEG-4');
    if links.total <= 60
        Movie.FrameRate = 5;
    elseif links.total > 60 && links.total <= 1000
        Movie.FrameRate = max([round(links.total/30); 10]);
    else
        Movie.FrameRate = min([round(links.total/60); 30]);
    end
    Movie.Quality = 100;
    open(Movie);
end

% Determines upper-limit value for the current calculations:
maxlim = max(max(abs(links.connections.currents)));
minlim = min(min(abs(nonzeros(links.connections.currents))));

% Define the color, line style, and line thickness:
gnd.color = [.75 .75 .75];

% Defining initial position of lightning link (m)
x1 = links.start(:,1)*sims.domain.dx;
y1 = links.start(:,2)*sims.domain.dy;
z1 = links.start(:,3)*sims.domain.dz+sims.domain.gnd; 

% Defining final position of lightning link (m)
x2 = links.end(:,1)*sims.domain.dx;
y2 = links.end(:,2)*sims.domain.dy;
z2 = links.end(:,3)*sims.domain.dz+sims.domain.gnd;

% Plots the cloud structure with the defined function below:
figure(1);
clf
set(gcf,'Position',[0,0,1.2*sims.plotWidth,1.2*sims.plotHeight]);
set(gcf,'Resize','off')

% Plotting link:
for ii = 1:1:links.total
    if strcmp(is.currentVisual,'log')
        intensity = reshape(((log10(abs(links.connections.currents(ii,1:ii))/minlim)))./(log10(maxlim/minlim)),[ii,1]);
    elseif strcmp(is.currentVisual,'linear')
        intensity = reshape(abs(links.connections.currents(ii,1:ii))./maxlim,[ii,1]);
    end
    coloring = ones([ii,3])-intensity;
    if links.connections.currents(ii,ii)<0
        coloring(1:ii,1) = 1; % red
    elseif links.connections.currents(ii,ii)>0
        coloring(1:ii,3) = 1; % blue
    end
    
    for N = 1:1:ii
        % Plotting initiation point:
        if N == 1
            if nodes.rho.values(ii,N)<0
                colorVector = [0 0 1];
            elseif nodes.rho.values(ii,N)>0
                colorVector = [1 0 0];
            elseif nodes.rho.values(ii,N)==0
                colorVector = [1 1 1];
            end
            % Ensures that a frame for the initiation point (step 0) is captured:
            if ii == 1
                scatter3(x1(N)*sims.spatialFactor.Number,y1(N)*sims.spatialFactor.Number,z1(N)*sims.spatialFactor.Number,50,'filled','Marker','*','MarkerEdgeColor',is.textColor,'MarkerFaceColor',is.textColor,'DisplayName','Initiation Point','LineWidth',2);
                hold on
                setUpAxes(sims,'xyz',is.backgroundColor);
                title(strcat("(", sims.objectType,") $n$ = ", int2str(ii-1) ," steps"),'FontSize',28,'FontWeight','bold','Interpreter','latex');
                legend('location','north','box','off','interpreter','latex','TextColor',is.textColor)
                if strcmp(is.Rec,'Y')
                    frame = getframe(gcf);
                    writeVideo(Movie,frame);
                    if ii == links.total
                        close(Movie);
                    end
                else
                    pause(0.01)
                end
                hold off
            end
            scatter3(x1(N)*sims.spatialFactor.Number,y1(N)*sims.spatialFactor.Number,z1(N)*sims.spatialFactor.Number,5+44*(abs(nodes.rho.intensity(ii,N))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
            hold on
            setUpAxes(sims,'xyz',is.backgroundColor);
        end
        % Determines the polarity of the charge density at the node:
        if nodes.rho.values(ii,N+1)<0
            colorVector = [0 0 1];
        elseif nodes.rho.values(ii,N+1)>0
            colorVector = [1 0 0];
        elseif nodes.rho.values(ii,N+1)==0
            colorVector = [1 1 1];
        end
        % Plots the links, nodes, and highlights the newest node:
        if N == 1 && ii == 1
            plot3([x1(N), x2(N)]*sims.spatialFactor.Number,[y1(N), y2(N)]*sims.spatialFactor.Number,[z1(N), z2(N)]*sims.spatialFactor.Number,'Color',coloring(:),'HandleVisibility','off','LineWidth',1);
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,0.1+49.9*(abs(nodes.rho.intensity(ii,N+1))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,5+49.9*(abs(nodes.rho.intensity(ii,N+1))),'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Current at tip: ",num2str(-links.connections.currents(N,N)*currentFactor.Number,'%+.2f')," ",currentFactor.LaTeX,"A"),'LineWidth',1);
            legend('location','north','box','off','interpreter','latex','TextColor',is.textColor)
        elseif N == ii            
            plot3([x1(N), x2(N)]*sims.spatialFactor.Number,[y1(N), y2(N)]*sims.spatialFactor.Number,[z1(N), z2(N)]*sims.spatialFactor.Number,'Color',coloring(N,:),'HandleVisibility','off','LineWidth',2);%,'LineWidth',1+4.*intensity(N));
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,0.1+49.9*(abs(nodes.rho.intensity(ii,N+1))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,5+49.9*(abs(nodes.rho.intensity(ii,N+1))),'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Current at tip: ",num2str(-links.connections.currents(N,N)*currentFactor.Number,'%+.2f')," ",currentFactor.LaTeX,"A"),'LineWidth',1);
            legend('location','north','box','off','interpreter','latex','TextColor',is.textColor)
        else
            plot3([x1(N), x2(N)]*sims.spatialFactor.Number,[y1(N), y2(N)]*sims.spatialFactor.Number,[z1(N), z2(N)]*sims.spatialFactor.Number,'Color',coloring(N,:),'HandleVisibility','off','LineWidth',1);%,'LineWidth',1+4.*intensity(N));
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,0.1+49.9*(abs(nodes.rho.intensity(ii,N+1))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
        end
    end
    % Sets the appropriate title:
    if ii == 1
        title(strcat("(", sims.objectType,") $n$ = ", int2str(ii) ," step $\rightarrow$ ",num2str(temporalFactor.Number*max(links.mintime(1:ii)),'%.2f')," ",temporalFactor.LaTeX,"s $\leq t_n \leq$ ",num2str(temporalFactor.Number*links.maxtime(ii),'%.2f')," ",temporalFactor.LaTeX,"s"),'FontSize',28,'FontWeight','bold','Interpreter','latex');
    else
        title(strcat("(", sims.objectType,") $n$ = ", int2str(ii) ," steps $\rightarrow$ ",num2str(temporalFactor.Number*max(links.mintime(1:ii)),'%.2f')," ",temporalFactor.LaTeX,"s $\leq t_n \leq$ ",num2str(temporalFactor.Number*links.maxtime(ii),'%.2f')," ",temporalFactor.LaTeX,"s"),'FontSize',28,'FontWeight','bold','Interpreter','latex');
    end
    % Save the frame if a video is requested:
    if strcmp(is.Rec,'Y')
        frame = getframe(gcf);
        writeVideo(Movie,frame);
        if ii == links.total
            close(Movie);
        end
    else
        pause(0.01)
    end
    % Exports the final figure:
    if ii == links.total
        if exist('export_fig') == 2 && strcmp(is.highResolution,'Y') == 1
            pause
            currentFolder = pwd;
            cd(sims.pathPNGs);
            export_fig HighRes_Current.png -transparent -m8
            cd(currentFolder);
        else
            exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentMovement_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
        end
    end
    hold off
end
