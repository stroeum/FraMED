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
clearvars -except sims is

%% Options to define variables that are typically input during run-time:
is.Rec             = 'N'; % Would you like to record the lightning propagation as a movie? (Y / N)
is.highResolution  = 'Y'; % Would you like to save the final image as a very high resolution image? (Y / N)
is.cropped         = 'Y'; % Would you like to automatically crop to the discharge? (Y / N)
is.plotRho         = 'N';
is.updateRho       = 'N';
is.Jefimenko       = 'Y';

is.plotWhich       = 'first'; % Provides the option to isolate plot by polarity. Options: 'all','pos','neg','first','final'
is.howMany         = 3;     % If <<is.plotWhich>> = 'first' || 'final', ensures only the <<is.plotWhich>> <<is.howMany>> frames will be plotted.

is.currentVisual   = 'log'; % Sets the scaling approach for the link's current intensity. Options: 'log' or 'linear'.
is.backgroundColor = 'w';   % Sets the plot's background color.
is.textColor       = 'k';   % Sets the plot's text color.

%% Fully-automated onwards:
fprintf('\n*** Executing Plot3D_TreeAndCurrents.m script. ***\n');  
if ~exist('sims','var')
    specifySimDetails;
end

% Read-in data for simulation and process the current at the tip:
[current, polarity, tau, links, rhos] = processTipCurrent(sims);

% Calculate the current through all links at all steps: 
[links, rhos] = processCurrents(links,rhos);

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
    fprintf(strcat("\t\t",num2str(sims.spatialFactor.Number*chargeTransported(end)./(links.pathdistance(end)/(polarity.neg(end)*sims.vprop.neg + polarity.pos(end)*sims.vprop.pos)),5),"\t",sims.spatialFactor.Unit,"A\t(initiation-to-ground path)\n"))
end

%% Plotting the current density:
% Defines the steps to be plotted:
switch is.plotWhich
    case 'all'
        loopArray = 1:1:links.total;
    case 'first'
        loopArray = 1:1:is.howMany;
    case 'final'
        loopArray = (links.total-is.howMany):1:links.total;
    case 'neg'
        loopArray = nonzeros((polarity.neg).*(1:1:links.total)')';
    case 'pos'
        loopArray = nonzeros((polarity.pos).*(1:1:links.total)')';
    otherwise
        fprintf('\nInvalid input for the is.plotWhich variable. Modify string and try again.\n');
        return
end
% Ensures initial and final states are plotted:
if loopArray(1)~=1 && ~strcmp(is.plotWhich,'final')
    loopArray = [1 loopArray];
end
if loopArray(end)~=links.total && ~strcmp(is.plotWhich,'first')
    loopArray = [loopArray links.total];
end

% Set movie recording
if (strcmp(is.Rec,'Y') == 1)
    if strcmp(is.plotWhich,'first') || strcmp(is.plotWhich,'final')
        Movie = VideoWriter(strcat(sims.pathVideos,'/TreeAndCurrents_',is.plotWhich,num2str(is.howMany),'_',sims.objectName,'_',sims.objectType,'_',is.cropped,is.plotRho),'MPEG-4');
    else
        Movie = VideoWriter(strcat(sims.pathVideos,'/TreeAndCurrents_',is.plotWhich,'_',sims.objectName,'_',sims.objectType,'_',is.cropped,is.plotRho),'MPEG-4');
    end
    totalFrames = length(loopArray);
    if totalFrames <= 60
        Movie.FrameRate = 2*(1+strcmp(is.plotWhich,'all'));
    elseif totalFrames > 60 && totalFrames <= 1000
        Movie.FrameRate = max([round(totalFrames/30); 2])*(1+strcmp(is.plotWhich,'all'));
    else
        Movie.FrameRate = min([round(totalFrames/60); 10])*(1+strcmp(is.plotWhich,'all'));
    end
    Movie.Quality = 100;
    open(Movie);
end
if strcmp(is.highResolution,'Y')
    sims.highRes = 'Y';
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

% Linear spaces for the three position dimensions:
x = (0:sims.domain.dx:sims.domain.maxx)'*sims.spatialFactor.Number;
y = (0:sims.domain.dy:sims.domain.maxy)'*sims.spatialFactor.Number;
z = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)'*sims.spatialFactor.Number;
[X,Y,Z]=meshgrid(x,y,z);
rho.data = convertTo3d(load("../results/rhoAmb.dat",'-ascii'),sims);
rho.max = max(max(max(rho.data)));
rho.min = min(min(min(rho.data)));

% Plots the cloud structure with the defined function below:
figure(1);
clf
set(gcf,'Position',[0,0,(1+strcmp(is.highResolution,'Y'))*(1+0.2*strcmp(is.plotRho,'Y'))*sims.plotWidth,(1+strcmp(is.highResolution,'Y'))*1.2*sims.plotHeight],'Color','White');
set(gcf,'Resize','off')


% Close-up Cropping:
if strcmp(is.cropped,'Y')
    ratio = max([max([max(links.start(:,1)) max(links.end(:,1))])./sims.domain.Nx; max([max(links.start(:,2)) max(links.end(:,2))])./sims.domain.Ny; max([max(links.start(:,3)) max(links.end(:,3))])./sims.domain.Nz;]);
else
    ratio = 1;
end
% Plotting link:
for ii = loopArray
    if polarity.pos(ii) == 1
        links.connections.currents(ii,1:ii) = -links.connections.currents(ii,1:ii);
        links.direction(ii,1:ii) = -links.direction(ii,1:ii);
    end
    polarity.step = ((links.connections.currents(ii,1:ii)>0)-(links.connections.currents(ii,1:ii)<0));
    links.connections.currents(ii,1:ii) = polarity.step.*links.connections.currents(ii,1:ii);
    links.direction(ii,1:ii) = polarity.step.*links.direction(ii,1:ii);

    if strcmp(is.currentVisual,'log')
        intensity = reshape(((log10(abs(links.connections.currents(ii,1:ii))/minlim)))./(log10(maxlim/minlim)),[ii,1]);
    elseif strcmp(is.currentVisual,'linear')
        intensity = reshape(abs(links.connections.currents(ii,1:ii))./maxlim,[ii,1]);
    end
    if strcmp(is.backgroundColor,'k')
        coloring = zeros([ii,3])+intensity;
    else
        coloring = ones([ii,3])-intensity;
    end
    % if ii==links.total && strcmp(sims.disType,'Cloud-to-Ground')
    %     coloring(links.connections.currents(ii,1:ii)<0,1) = 1; % red
    %     coloring(links.connections.currents(ii,1:ii)>0,3) = 1; % blue
    % end
    for N = 1:1:ii
        % Plotting initiation point:
        if N == 1
            if rhos.values(ii,N)<0
                colorVector = [0 0 1];
            elseif rhos.values(ii,N)>0
                colorVector = [1 0 0];
            elseif rhos.values(ii,N)==0
                colorVector = [1 1 1];
            end
            % Ensures that a frame for the initiation point (step 0) is captured:
            if ii == 1
                scatter3(0,0,0,(1/ratio)*(1+3*strcmp(is.highResolution,'Y'))*50,'filled','Marker','square','MarkerEdgeColor',is.textColor,'MarkerFaceColor',is.textColor,'DisplayName','Ground Station');
                hold on
                if strcmp(is.plotRho,'Y')
                    plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);
                end
                if strcmp(is.highResolution,'Y')
                    setUpAxes(sims,'xyz_hires',is.backgroundColor);
                else
                    setUpAxes(sims,'xyz',is.backgroundColor);
                end
                
                if strcmp(is.cropped,'Y')
                    xlim([0 ceil(ratio*sims.domain.dx*sims.domain.Nx*sims.spatialFactor.Number)]);
                    ylim([0 ceil(ratio*sims.domain.dy*sims.domain.Ny*sims.spatialFactor.Number)]);
                    zlim([0 ceil(ratio*sims.domain.dz*sims.domain.Nz*sims.spatialFactor.Number)]);
                end
                scatter3(x1(N)*sims.spatialFactor.Number,y1(N)*sims.spatialFactor.Number,z1(N)*sims.spatialFactor.Number,(1/ratio)*(1+3*strcmp(is.highResolution,'Y'))*35,'filled','Marker','*','MarkerEdgeColor',is.textColor,'MarkerFaceColor',is.textColor,'DisplayName','Initiation Point','LineWidth',(1+strcmp(is.highResolution,'Y'))*1.5);
                sgtitle(strcat("\textbf{Ambient State}"),'FontSize',(1+strcmp(is.highResolution,'Y'))*30,'FontWeight','bold','Interpreter','latex');
                title({" "},{strcat("($n$ = ", int2str(ii-1) ," steps)")},'FontSize',(1+strcmp(is.highResolution,'Y'))*26,'FontWeight','bold','Interpreter','latex');
                [~,objh] = legend('location','north','box','off','interpreter','latex','TextColor',is.textColor,'FontSize',(1+strcmp(is.highResolution,'Y'))*16);
                markerMod = findobj(objh,'Type','patch');
                set(markerMod(1),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*10);
                set(markerMod(2),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*6);
                if strcmp(is.Rec,'Y')
                    frame = getframe(gcf);
                    writeVideo(Movie,frame);
                    if ii == loopArray(end)
                        close(Movie);
                    end
                else
                    pause(0.01)
                end
                exportgraphics(gcf,strcat(sims.pathPNGs,'/InitialState_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
                sgtitle(strcat("\textbf{",sims.disType," ",sims.objectType," on ",sims.objectName,"}"),'FontSize',(1+strcmp(is.highResolution,'Y'))*30,'FontWeight','bold','Interpreter','latex');
                hold off
            end
            scatter3(0,0,0,(1/ratio)*(1+3*strcmp(is.highResolution,'Y'))*50,'filled','Marker','square','MarkerEdgeColor',is.textColor,'MarkerFaceColor',is.textColor,'DisplayName','Ground Station');
            hold on
            scatter3(x1(N)*sims.spatialFactor.Number,y1(N)*sims.spatialFactor.Number,z1(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(5+44*(abs(rhos.intensity(ii,N)))),'filled','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'HandleVisibility','off');
        
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
                    if strcmp(is.plotRho,'Y')
                        plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);
                    end
                    if strcmp(is.highResolution,'Y')
                        setUpAxes(sims,'xyz_hires',is.backgroundColor);
                    else
                        setUpAxes(sims,'xyz',is.backgroundColor);
                    end
                end
            else
                if strcmp(is.plotRho,'Y')
                    plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);
                else
                    if strcmp(is.highResolution,'Y')
                        setUpAxes(sims,'xyz_hires',is.backgroundColor);
                    else
                        setUpAxes(sims,'xyz',is.backgroundColor);
                    end
                end
            end

            if strcmp(is.cropped,'Y')
                xlim([0 ceil(ratio*sims.domain.dx*sims.domain.Nx*sims.spatialFactor.Number)]);
                ylim([0 ceil(ratio*sims.domain.dy*sims.domain.Ny*sims.spatialFactor.Number)]);
                zlim([0 ceil(ratio*sims.domain.dz*sims.domain.Nz*sims.spatialFactor.Number)]);
            end
        end
        % Determines the polarity of the charge density at the node:
        if rhos.values(ii,N+1)<0 || (ii==N && polarity.neg(ii)==1)
            colorVector = [0 0 1];
        elseif rhos.values(ii,N+1)>0 || (ii==N && polarity.pos(ii)==1)
            colorVector = [1 0 0];
        elseif rhos.values(ii,N+1)==0
            colorVector = [1 1 1];
        end
        % Plots the links, nodes, and highlights the newest node:
        if N == 1 && ii == 1
            if links.direction(ii,N) == 1
                mArrow3([x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,[x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,'color',coloring(:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            else
                mArrow3([x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,[x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,'color',coloring(:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            end
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(0.1+49.9*(abs(rhos.intensity(ii,N+1)))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(5+49.9*(abs(rhos.intensity(ii,N+1)))),'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Current at tip: ",num2str(abs(links.connections.currents(N,N))*currentFactor.Number,'%.3f')," ",currentFactor.LaTeX,"A"),'LineWidth',0.01+(abs(rhos.intensity(ii,N+1))));
            [~,objh] = legend('location','north','box','off','interpreter','latex','TextColor',is.textColor,'FontSize',(1+strcmp(is.highResolution,'Y'))*16);
            markerMod = findobj(objh,'Type','patch');
            set(markerMod(1),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*10);
            set(markerMod(2),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*6);
        elseif N == ii            
            if intensity(ii) == 0
                coloring(N,:) = [0 0 0];
                intensity(N) = 1;
            end
            if links.direction(ii,N) == 1
                mArrow3([x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,[x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,'color',coloring(N,:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            else
                mArrow3([x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,[x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,'color',coloring(N,:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            end
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(0.1+49.9*(abs(rhos.intensity(ii,N+1)))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
            if ii==links.total && (strcmp(sims.disType,'Horizontal') || strcmp(sims.disType,'Jet'))
                scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(5+49.9*(abs(rhos.intensity(ii,N+1)))),'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Interacts with Boundary"),'LineWidth',(1+strcmp(is.highResolution,'Y'))*(0.01+(abs(rhos.intensity(ii,N+1)))));
            elseif ii==links.total && strcmp(sims.disType,'Cloud-to-Ground') 
                scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(1/ratio)*(1+3*strcmp(is.highResolution,'Y'))*25,'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Current to ground: ",num2str(abs(links.connections.currents(N,N))*currentFactor.Number,'%.3f')," ",currentFactor.LaTeX,"A"),'LineWidth',(1+strcmp(is.highResolution,'Y'))*(0.01+(abs(rhos.intensity(ii,N+1)))));
            else
                scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(5+49.9*(abs(rhos.intensity(ii,N+1)))),'Marker','o','MarkerEdgeColor',is.textColor,'MarkerFaceColor',colorVector,'DisplayName',strcat("Current at tip: ",num2str(abs(links.connections.currents(N,N))*currentFactor.Number,'%.3f')," ",currentFactor.LaTeX,"A"),'LineWidth',(1+strcmp(is.highResolution,'Y'))*(0.01+(abs(rhos.intensity(ii,N+1)))));
            end
            [~,objh] = legend('location','north','box','off','interpreter','latex','TextColor',is.textColor,'FontSize',(1+strcmp(is.highResolution,'Y'))*16);
            markerMod = findobj(objh,'Type','patch');
            set(markerMod(1),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*10);
            set(markerMod(2),'MarkerSize',(1+strcmp(is.highResolution,'Y'))*6);
        else
            if links.direction(ii,N) == 1
                mArrow3([x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,[x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,'color',coloring(N,:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            else
                mArrow3([x2(N) y2(N) z2(N)].*sims.spatialFactor.Number,[x1(N) y1(N) z1(N)].*sims.spatialFactor.Number,'color',coloring(N,:),'stemWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.03*intensity(N),'tipWidth',(1+0.01*strcmp(is.cropped,'Y'))*0.05,'HandleVisibility','off');
            end
            scatter3(x2(N)*sims.spatialFactor.Number,y2(N)*sims.spatialFactor.Number,z2(N)*sims.spatialFactor.Number,(0.5/ratio)*(1+3*strcmp(is.highResolution,'Y'))*(0.1+49.9*(abs(rhos.intensity(ii,N+1)))),'filled','MarkerEdgeColor',colorVector,'MarkerFaceColor',colorVector,'HandleVisibility','off');
        end
    end
    
    % Reverts values to there original state:
    links.connections.currents(ii,1:ii) = polarity.step.*links.connections.currents(ii,1:ii);
    links.direction(ii,1:ii) = polarity.step.*links.direction(ii,1:ii);
    if polarity.pos(ii) == 1
        links.connections.currents(ii,1:ii) = -links.connections.currents(ii,1:ii);
        links.direction(ii,1:ii) = -links.direction(ii,1:ii);
    end

    % Sets the appropriate title:
    if ii == 1
        title({" "},{strcat("($n$ = ", int2str(ii) ," step $\rightarrow$ ",num2str(temporalFactor.Number*max(links.mintime(1:ii)),'%.3f')," ",temporalFactor.LaTeX,"s $\leq t_n \leq$ ",num2str(temporalFactor.Number*links.maxtime(ii),'%.3f')," ",temporalFactor.LaTeX,"s)")},'FontSize',(1+strcmp(is.highResolution,'Y'))*26,'FontWeight','bold','Interpreter','latex');
    else
        if ii == links.total
            sgtitle(strcat("\textbf{Final State}"),'FontSize',(1+strcmp(is.highResolution,'Y'))*30,'FontWeight','bold','Interpreter','latex');
        elseif strcmp(is.plotWhich,'final') && ii == loopArray(1)
            sgtitle(strcat("\textbf{",sims.disType," ",sims.objectType," on ",sims.objectName,"}"),'FontSize',(1+strcmp(is.highResolution,'Y'))*30,'FontWeight','bold','Interpreter','latex');
        end
        title({" "},{strcat("($n$ = ", int2str(ii) ," steps $\rightarrow$ ",num2str(temporalFactor.Number*max(links.mintime(1:ii)),'%.3f')," ",temporalFactor.LaTeX,"s $\leq t_n \leq$ ",num2str(temporalFactor.Number*links.maxtime(ii),'%.3f')," ",temporalFactor.LaTeX,"s)")},'FontSize',(1+strcmp(is.highResolution,'Y'))*26,'FontWeight','bold','Interpreter','latex');
    end
    
    % Save the frame if a video is requested:
    if strcmp(is.Rec,'Y')
        frame = getframe(gcf);
        writeVideo(Movie,frame);
        if ii == loopArray(end)
            close(Movie);
        end
    else
        exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentMovement_Step-',int2str(ii),'_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
                
        pause(0.01)
    end
    % Exports the final figure:
    if ii == loopArray(end)
        exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentMovement_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
        % if exist('export_fig') == 2 && strcmp(is.highResolution,'Y') == 1
        %     pause
        %     currentFolder = pwd;
        %     cd(sims.pathPNGs);
        %     export_fig HighRes_Current.png -transparent -m8
        %     cd(currentFolder);
        % else
        %     exportgraphics(gcf,strcat(sims.pathPNGs,'/CurrentMovement_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',600);
        % end
    end
    hold off
end

if strcmp(is.highResolution,'Y')
    sims.highRes = 'N';
end

%% Georg Stillfried (2025). mArrow3.m - easy-to-use 3D arrow 
% (https://www.mathworks.com/matlabcentral/fileexchange/25372-marrow3-m-easy-to-use-3d-arrow), 
% MATLAB Central File Exchange. Retrieved December 10, 2025.
function h = mArrow3(p1,p2,varargin)
%mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
%
% syntax:   h = mArrow3(p1,p2)
%           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
%
% with:     p1:         starting point
%           p2:         end point
%           properties: 'color':      color according to MATLAB specification
%                                     (see MATLAB help item 'ColorSpec')
%                       'stemWidth':  width of the line
%                       'tipWidth':   width of the cone                       
%
%           Additionally, you can specify any patch object properties. (For
%           example, you can make the arrow semitransparent by using
%           'facealpha'.)
%                       
% example1: h = mArrow3([0 0 0],[1 1 1])
%           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
%
% example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
%           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
%
% hint:     use light to achieve 3D impression
%
    propertyNames = {'edgeColor'};
    propertyValues = {'none'};    
    % evaluate property specifications
    for argno = 1:2:nargin-2
        switch varargin{argno}
            case 'color'
                propertyNames = {propertyNames{:},'facecolor'};
                propertyValues = {propertyValues{:},varargin{argno+1}};
            case 'stemWidth'
                if isreal(varargin{argno+1})
                    stemWidth = varargin{argno+1};
                else
                    warning('mArrow3:stemWidth','stemWidth must be a real number');
                end
            case 'tipWidth'
                if isreal(varargin{argno+1})
                    tipWidth = varargin{argno+1};
                else
                    warning('mArrow3:tipWidth','tipWidth must be a real number');
                end
            otherwise
                propertyNames = {propertyNames{:},varargin{argno}};
                propertyValues = {propertyValues{:},varargin{argno+1}};
        end
    end            
    % default parameters
    if ~exist('stemWidth','var')
        ax = axis;
        if numel(ax)==4
            stemWidth = norm(ax([2 4])-ax([1 3]))/300;
        elseif numel(ax)==6
            stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
        end
    end
    if ~exist('tipWidth','var')
        tipWidth = 3*stemWidth;
    end
    tipAngle = 22.5/180*pi;
    tipLength = tipWidth/tan(tipAngle/2);
    ppsc = 50;  % (points per small circle)
    ppbc = 250; % (points per big circle)
    % ensure column vectors
    p1 = p1(:);
    p2 = p2(:);
    % basic lengths and vectors
    x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
    y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
    if norm(y)<0.1
        y = cross(x,[0;1;0]);
    end
    y = y/norm(y);
    z = cross(x,y);
    z = z/norm(z);
    % basic angles
    theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
    sintheta = sin(theta);
    costheta = cos(theta);
    upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
    sinupsilon = sin(upsilon);
    cosupsilon = cos(upsilon);
    % initialize face matrix
    f = NaN([ppsc+ppbc+2 ppbc+1]);
    % normal arrow
    if norm(p2-p1)>tipLength
        % vertices of the first stem circle
        for idx = 1:ppsc+1
            v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
        end
        % vertices of the second stem circle
        p3 = p2-tipLength*x;
        for idx = 1:ppsc+1
            v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
        end
        % vertices of the tip circle
        for idx = 1:ppbc+1
            v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
        end
        % vertex of the tiptip
        v(2*ppsc+ppbc+4,:) = p2;
        % face of the stem circle
        f(1,1:ppsc+1) = 1:ppsc+1;
        % faces of the stem cylinder
        for idx = 1:ppsc
            f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
        end
        % face of the tip circle
        f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
        % faces of the tip cone
        for idx = 1:ppbc
            f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
        end
    % only cone v
    else
        tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
        % vertices of the tip circle
        for idx = 1:ppbc+1
            v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
        end
        % vertex of the tiptip
        v(ppbc+2,:) = p2;
        % face of the tip circle
        f(1,:) = 1:ppbc+1;
        % faces of the tip cone
        for idx = 1:ppbc
            f(1+idx,1:3) = [idx idx+1 ppbc+2];
        end
    end
    % draw
    fv.faces = f;
    fv.vertices = v;
    h = patch(fv);
    for propno = 1:numel(propertyNames)
        try
            set(h,propertyNames{propno},propertyValues{propno});
        catch
            disp(lasterr)
        end
    end
end