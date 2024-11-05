% Initiate
close all
clearvars -except sims
clf

% For formatting exported image size:
vals.finalfieldflag = 0;

%% Provides the option to define variables that are typically input during run-time:
%is.Rec            = 'Y'; % Would you like to record the lightning propagation as a movie? (Y / N)
%is.updateRho      = 'Y'; % Would you like to update the charge distribution coloring for every saved step? (Y / N)
%is.highResolution = 'N'; % Would you like to save the final image as a very high resolution image? (Y / N)

if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
    prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
    sims.objectType = input(prompt2,'s');
    while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
        fprintf('\tNot an acceptable input. Please enter Streamer or Leader.\n');
        sims.objectType = input(prompt2,'s');
    end

    % Settings to ensure proper directory referencing:
    sims.pathPNGs = ['../Figures/',sims.objectName,'/',sims.objectType,'/PNGs'];
    if ~exist(sims.pathPNGs,'dir')
        mkdir(sims.pathPNGs);
    end
    sims.pathVideos = ['../Figures/',sims.objectName,'/',sims.objectType,'/Videos'];
    if ~exist(sims.pathVideos,'dir')
        mkdir(sims.pathVideos);
    end

    % Specifies the boundary conditions for the simulation:
    prompt_BCtype = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
    sims.BCtype = input(prompt_BCtype,'s');                    
    while ~strcmp(sims.BCtype,'FS') && ~strcmp(sims.BCtype,'G')
        fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
        sims.BCtype = input(prompt_BCtype,'s');
    end
end 

cd ../results

% Load data files
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
load('InitPoint.dat',        '-ascii');
rho.data  = load('rhoAmb.dat',           '-ascii');
Links.ID  = load('EstablishedLinks.dat', '-ascii');
gnd.alt   = load('z_gnd.dat',            '-ascii');
stepsaves = abs(load('step3d.dat',       '-ascii'));
ChargeLayers = load('ChargeLayers.dat');

if isempty(Links.ID)
    fprintf('\n*** LightningVisual.m cannot be executed with current EstablishedLinks.dat file. ***\n');
    cd ../viz
    return
elseif ~exist('is','var')
    fprintf('\n*** Executing LightningVisual.m script. ***\n');
    % User-Based Settings:
    if ~isfield(is,'Rec')
        prompt_Rec = '\nWould you like to record the lightning propagation as a movie? (Y / N)\n-->';
        is.Rec = input(prompt_Rec,'s');                    
        while ~strcmp(is.Rec,'Y') && ~strcmp(is.Rec,'N')
            fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
            is.Rec = input(prompt_Rec,'s');
        end
    end
    if ~isfield(is,'updateRho')
        prompt_updateRho = '\nWould you like to update the charge distribution coloring for every saved step? (Y / N)\n-->';
        is.updateRho = input(prompt_updateRho,'s');                    
        while ~strcmp(is.updateRho,'Y') && ~strcmp(is.updateRho,'N')
            fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
            is.updateRho = input(prompt_updateRho,'s');
        end
    end
    if ~isfield(is,'highResolution')
        prompt_highResolution = '\nWould you like to save the final image as a very high resolution image? (Y / N)\nWARNING: Only recommended for preparing posters.\n-->';
        is.highResolution = input(prompt_highResolution,'s');                    
        while ~strcmp(is.highResolution,'Y') && ~strcmp(is.highResolution,'N')
            fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
            is.highResolution = input(prompt_highResolution,'s');
        end
    end
end

% Derive main parameters
N.x = Nxyz(1);
N.y = Nxyz(2);
N.z = Nxyz(3);

d.x = dxyz(1);           % _m
d.y = dxyz(2);           % _m
d.z = dxyz(3);           % _m

L.x = (N.x-1)*d.x;         % _m
L.y = (N.y-1)*d.y;         % _m
L.z = (N.z-1)*d.z;         % _m

S.x = InitPoint(1);   % _m
S.y = InitPoint(2);   % _m
S.z = InitPoint(3);   % _m
S.R = InitPoint(4);   % _m

S.i = round(S.x/d.x);
S.j = round(S.y/d.y);
S.k = round(S.z/d.z);

dx = dxyz(1);               % _m
dy = dxyz(2);               % _m
dz = dxyz(3);               % _m

Lx = (Nxyz(1)-1)*dx*1e-3;         % _km
Ly = (Nxyz(2)-1)*dy*1e-3;         % _km
Lz = (Nxyz(3)-1)*dz*1e-3;         % _km

%% Derive data for plotting
x = ((0:(N.x-1))*d.x)*1e-3;
y = ((0:(N.y-1))*d.y)*1e-3;
z = ((0:(N.z-1))*d.z + gnd.alt)*1e-3;
[vals.y,vals.z1]    = meshgrid(y,z);
[vals.x,vals.z2]    = meshgrid(x,z);
Links.Nb = size(Links.ID);
Links.Nb = Links.Nb(1);

% Assigns plot height and width based on domain:
if ~isfield(sims,'plotWidth') || ~isfield(sims,'plotHeight')
    sims.plotWidth = 600;
    sims.plotHeight = round((Nxyz(3)/max(Nxyz(1:2)))*5)*80;
end

clear dxyz
clear InitPoint
cd ../viz

rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3

[X,Y,Z] = meshgrid(x,y,z);
rho.max = .95* max(max(max(rho.data)));
rho.min = .95* min(min(min(rho.data)));

% Map ColorScale
gnd.color = [.75 .75 .75];
color     = colormap(jet(Links.Nb));
is.BW = 1; %input('Is plot Monochrome? (1: yes, else: no)\n>> ');
if (is.BW == 1)
    for ii=1:Links.Nb
        color(ii,:) = [0 0 0];
    end
end
% Set movie recording
if (strcmp(is.Rec,'Y') == 1)
    Movie = VideoWriter([sims.pathVideos,'/',sims.objectName,'_',sims.objectType,'FieldLinksVideo'],'MPEG-4');
    if Links.Nb <= 60
        Movie.FrameRate = 1;
    elseif Links.Nb > 60 && Links.Nb <= 1000
        Movie.FrameRate = round(Links.Nb/30);
    else
        Movie.FrameRate = round(Links.Nb/60);
    end
    Movie.Quality = 100;
    open(Movie);
end


% Draw the tree
subplot(1,3,2)
hold on;
set(gcf,'Position',[0,0,sims.plotWidth*2,sims.plotHeight/1.5]);
set(gcf,'Resize','off');
grid on;

% Plots the cloud structure with the defined function below:
axis equal
xlim([0 max(x)]);
ylim([0 max(y)]);
zlim([0 max(z)]);

% Represents the neutrally charged (grounded) surface:
if strcmp(sims.BCtype,'G')
    P.x = [L.x 0 0 L.x]*1e-3;
    P.y = [L.y L.y 0 0]*1e-3;
    P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
    patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
end
hold off
        
% Initialize distance traveled for lightning links:
distance = 0;
% For-loop to plot lightning discharge links:
for ii=1:Links.Nb
    if mod((ii-1),stepsaves) == 0
        % Determines the 2D electric field values for the step:
        E = consolidateEfield((ii-1),vals,Nxyz);
        if ii == 1
            % Plots the charge density distribution in panel 2:
            subplot(1,3,2)
            hold on
            chargeDensityUpdate('white',0.25,rho,X,Y,Z);
            ax2 = gca;
            hold off
                            
            % Plots the electric field lines and ensures proper formatting for panel 3:
            subplot(1,3,3)
            hold on
            subplot3Position = get(subplot(1,3,3),'Position');
            sp3.x = subplot3Position(1);
            sp3.y = subplot3Position(2);
            sp3.dx = subplot3Position(3);
            sp3.dy = subplot3Position(4);
            p3 = streamslice(vals.x,vals.z2,E.x2D,E.z2Dx,10,'arrows');
            set(findobj('Type','line'),'Color','k')
            sub3Inner = get(subplot(1,3,3),'InnerPosition');
            set(subplot(1,3,3),'Position',sub3Inner,'PositionConstraint','innerposition');
            set(gca,'FontSize',10,'TickLabelInterpreter','latex','YAxisLocation','right');
            xlabel('$x$ (km)','Interpreter','latex','FontSize',16);
            ylabel('$z$ (km)','Interpreter','latex','FontSize',16);
            axis xy
            axis image
            set(gca,'XLim',[min(vals.x(:)) max(vals.x(:))],'YLim',[min(vals.z2(:)) max(vals.z2(:))]);
            ax3 = gca;
            box on
            hold off
            sub2Pos = get(subplot(1,3,2),'Position');
            sub2Border = get(subplot(1,3,2),'OuterPosition');
            sub3Outer = get(subplot(1,3,3),'OuterPosition');

            % Plots the electric field lines and ensures proper formatting for panel 1:
            subplot(1,3,1)
            hold on
            subplot1Position = get(subplot(1,3,1),'Position');
            sp1.x = subplot1Position(1);
            sp1.y = subplot1Position(2);
            sp1.dx = subplot1Position(3);
            sp1.dy = subplot1Position(4);
            p1 = streamslice(vals.y,vals.z1,E.y2D,E.z2Dy,10,'arrows');
            set(findobj('Type','line'),'Color','k');
            sub1Inner = [(sub2Border(1)-(sub3Inner(1)-(sub2Border(1)+sub2Border(3)))-(sub3Inner(3)*(Nxyz(2)/Nxyz(1)))) sub3Inner(2) (sub3Inner(3)*(Nxyz(2)/Nxyz(1))) sub3Inner(4)];
            set(subplot(1,3,1),'Position',sub1Inner,'PositionConstraint','innerposition');
            set(gca,'FontSize',10,'TickLabelInterpreter','latex');
            xlabel('$y$ (km)','Interpreter','latex','FontSize',16);
            ylabel('$z$ (km)','Interpreter','latex','FontSize',16);
            axis xy
            axis image
            set(gca,'XDir','reverse','XLim',[min(vals.y(:)) max(vals.y(:))],'YLim',[min(vals.z1(:)) max(vals.z1(:))]);
            ax1 = gca;
            box on
            hold off
        else % If it is a multiple of the step save extracts, just not the initial:
            hold off

            % Updates the field lines and reformats panel 1:
            subplot(1,3,1)
            hold on
            delete(p1);
            p1 = streamslice(vals.y,vals.z1,E.y2D,E.z2Dy,10,'arrows');
            set(findobj('Type','line'),'Color','k');
            set(subplot(1,3,1),'Position',sub1Inner,'PositionConstraint','innerposition');
            set(gca,'FontSize',10,'TickLabelInterpreter','latex');
            xlabel('$y$ (km)','Interpreter','latex','FontSize',16);
            ylabel('$z$ (km)','Interpreter','latex','FontSize',16);
            axis xy
            axis image
            set(gca,'XDir','reverse','XLim',[min(vals.y(:)) max(vals.y(:))],'YLim',[min(vals.z1(:)) max(vals.z1(:))]);
            ax1 = gca;
            box on
            hold off

            % Updates the field lines and reformats panel 3:
            subplot(1,3,3)
            delete(p3);
            hold on
            p3 = streamslice(vals.x,vals.z2,E.x2D,E.z2Dx,10,'arrows');
            set(findobj('Type','line'),'Color','k');
            set(subplot(1,3,3),'Position',sub3Inner,'PositionConstraint','innerposition');
            set(gca,'FontSize',10,'TickLabelInterpreter','latex','YAxisLocation','right');
            xlabel('$x$ (km)','Interpreter','latex','FontSize',16);
            ylabel('$z$ (km)','Interpreter','latex','FontSize',16);
            axis xy
            axis image
            set(gca,'XLim',[min(vals.x(:)) max(vals.x(:))],'YLim',[min(vals.z2(:)) max(vals.z2(:))]);
            ax3 = gca;
            box on
            hold off

            % Updates the charge density distribution in panel 2
            subplot(1,3,2)
            hold on
            allPatches = findall(subplot(1,3,2),'type','patch');
            delete(allPatches);
            if strcmp(sims.BCtype,'G')
                P.x = [L.x 0 0 L.x]*1e-3;
                P.y = [L.y L.y 0 0]*1e-3;
                P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
                patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
            end
            rho.data = load(['../results/rho3d',num2str(ii-1),'.dat'],           '-ascii');
            rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3
            rho.max = .95* max(max(max(rho.data)));
            rho.min = .95* min(min(min(rho.data)));
            chargeDensityUpdate('white',0.4,rho,X,Y,Z);
            hold off
        end

    end  
    % Revisiting panel 2 to draw the link propagations:
    subplot(1,3,2)
    hold on

    % Defining initial position of lightning link (m)
    x1 = Links.ID(ii,1)*d.x;
    y1 = Links.ID(ii,2)*d.y;
    z1 = Links.ID(ii,3)*d.z+gnd.alt; 

    % Defining final position of lightning link (m)
    x2 = Links.ID(ii,4)*d.x;
    y2 = Links.ID(ii,5)*d.y;
    z2 = Links.ID(ii,6)*d.z+gnd.alt;

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

    % Plotting link:
    plot3(...
        [x1, x2]*1e-3,...
        [y1, y2]*1e-3,...
        [z1, z2]*1e-3,...
        'Color',color(ii,:),'HandleVisibility','off');
    
    % Formatting title to reflect step value:
    if ii == 1
        titleformat = sgtitle(['(',sims.objectName,'): ',sims.objectType,' discharge after ', int2str(ii) ,' step(s)'],'FontSize',24,'FontWeight','bold','Interpreter','latex');
        set(gcf,'Position',[0,0,sims.plotWidth*2,sims.plotHeight/1.5]);
        set(gcf,'Resize','off');
    else
        delete(titleformat);
        titleformat = sgtitle(['(',sims.objectName,'): ',sims.objectType,' discharge after ', int2str(ii) ,' step(s)'],'FontSize',24,'FontWeight','bold','Interpreter','latex');
        set(gcf,'Position',[0,0,sims.plotWidth*2,sims.plotHeight/1.5]);
        set(gcf,'Resize','off');
    end
    box on
    if mod((ii-1),stepsaves) == 0
        exportgraphics(gcf,[sims.pathPNGs,'/LightningFields_',num2str(ii-1),'_',sims.objectName,'_',sims.objectType,'.png'],'Resolution',600);
    end
    if(strcmp(is.Rec,'Y') == 1)
        set(gcf,'Position',[0,0,sims.plotWidth*2,sims.plotHeight/1.5]);
        set(gcf,'Resize','off');
        frame = getframe(gcf);
        writeVideo(Movie,frame);
    end
end
fprintf(['\n\t',sims.objectType,' reaches minimum of ',num2str(minHeight,'%.2f'),' meters (',num2str(minHeight-initHeight,'%+.2f'),' meters)\n']);
fprintf(['\n\t',sims.objectType,' initiated at ',num2str(initHeight,'%.2f'),' meters\n']);
fprintf(['\n\t',sims.objectType,' reaches maximum of ',num2str(maxHeight,'%.2f'),' meters (',num2str(maxHeight-initHeight,'%+.2f'),' meters)\n']);
fprintf(['\n\t',sims.objectType,' has propagated a total of ',num2str(distance,'%.2f'),' meters\n']);
hold off;
% Record the movie
if (strcmp(is.Rec,'Y') == 1)
    set(gcf,'Position',[0,0,sims.plotWidth*2,sims.plotHeight/1.5]);
    set(gcf,'Resize','off');
    frame = getframe(gcf);
    writeVideo(Movie,frame);
    close(Movie);
end

% Exports the final graphic:
exportgraphics(gcf,[sims.pathPNGs,'/LightningFields_',sims.objectName,'_',sims.objectType,'.png'],'Resolution',600);
    
%% Functions:
% Converts 3D electric field data files into respective 2D components:
function E = consolidateEfield(num,vals,Nxyz)
    if vals.finalfieldflag == 1
        Ex3D         = ConvertTo3d(load('../results/Ex3d.dat'),Nxyz);
        Ey3D         = ConvertTo3d(load('../results/Ey3d.dat'),Nxyz);
        Ez3D         = ConvertTo3d(load('../results/Ez3d.dat'),Nxyz);
    else
        Ex3D         = ConvertTo3d(load(['../results/Ex3d',num2str(num),'.dat']),Nxyz);
        Ey3D         = ConvertTo3d(load(['../results/Ey3d',num2str(num),'.dat']),Nxyz);
        Ez3D         = ConvertTo3d(load(['../results/Ez3d',num2str(num),'.dat']),Nxyz);
    end
    E.x2D          = permute(Ex3D(:,(Nxyz(2)-1)/2,:),[3,1,2]);
    E.y2D          = permute(Ey3D((Nxyz(1)-1)/2,:,:),[3,2,1]);
    E.z2Dx         = permute(Ez3D(:,(Nxyz(2)-1)/2,:),[3,1,2]);
    E.z2Dy         = permute(Ez3D((Nxyz(1)-1)/2,:,:),[3,2,1]);
    clear Ex3D Ey3D Ez3D
end

% Similar functionality to plottingChargeRegions.m file, but with custom
% colorbar placement for a compact form:
function chargeDensityUpdate(colorbarRange,alphaValue,rhoDataOG,Xval,Yval,Zval)
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
    if size(findall(gcf,'type','colorbar'),1) == 0
        c = colorbar;
        c.Label.String = 'Charge Density (nC/m$^3$)';
        c.Label.Interpreter = 'latex';
        c.Label.HorizontalAlignment = 'center';
        c.Label.FontSize = 16;
        c.TickDirection = 'out';
        c.TickLabelInterpreter = 'latex';
        c.AxisLocation = "out";
        c.Location = "southoutside"; 
    end
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 10;
    xlabel('$x$ (km)','Interpreter','latex','FontSize',10,'HorizontalAlignment','left');
    ylabel('$y$ (km)','Interpreter','latex','FontSize',10,'HorizontalAlignment','right');
    grid on
    view(-45,9)
end