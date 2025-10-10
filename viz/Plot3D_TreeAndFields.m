% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: Plot3D_TreeAndFields.m                                      %
%    Purpose: Creates a three panel graphic that plots the electric field %
%             lines for requested perspectives next to the propagating    %
%             discharge. Discharge can be overlaid with a visualization   %
%             of the potential or charge density at the associated step.  %
%             Can be turned into a video as well.                         %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: March 11, 2024                                              %
%    Updates: February 2025 - Introduced a custom colormap for field-like %
%                             values, the ability to zoom-in, XY plane    %
%                             extracts, and the choice between potential  %
%                             and charge density overlays.                % 
%              October 2025 - Integrated the checkMagnitude.m function    %
%                             for plot formatting and axes labeling.      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clearvars -except sims

%% User-defined variables that influence the plot:
is.Rec            = 'N'; % Would you like to record the lightning propagation as a video? (Y / N)
is.monoChrome     = 'N'; % Would you like the discharge plotted only in black? (Y / N)
is.whichScalar    = 'P'; % Would you like to plot the charge density (C) or the potential (P)?
is.updateScalar   = 'N'; % Would you like to update the charge/potential distribution coloring using the saved steps? (Y / N)
is.highRes        = 'N';

is.zoomedIn       = 'N'; % Would you like to plot a specific region of the domain (Y/N)?
is.zoomedX        = [7500; 17500];  % (is.zoomedIn == 'Y'): Zoomed in region for the X domain, in meters.
is.zoomedY        = [7500; 17500];  % (is.zoomedIn == 'Y'): Zoomed in region for the Y domain, in meters.
is.zoomedZ        = [4000; 14000];  % (is.zoomedIn == 'Y'): Zoomed in region for the Z domain, in meters.

is.horizontal     = 'N'; % Would you like to plot the field for the XY plane at two altitudes (Y) or for the XZ/YZ planes at the center of the domain (N)?
is.altitudes      = [5000; 6000];  % (is.horizontal == 'Y'): Altitudes for horizontal plane plots, in meters.

is.singleImages   = 'N'; % Would you like to export single images for all of the available saved steps? (Y / N)
is.setExtracts    = []; % Would you like to export for any explicit steps? 

%% Fully-automated onwards:
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    specifySimDetails;
end

% Load data files
cd ../results
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
load('MaximumEfield.dat',    '-ascii');
Links.ID  = load('EstablishedLinks.dat', '-ascii');
gnd.alt   = load('z_gnd.dat',            '-ascii');
stepsaves = abs(load('step3d.dat',       '-ascii'));
polarity  = load('TransportedRhoEnd.dat','-ascii');
cd ../viz

if isempty(Links.ID)
    fprintf('\n*** Plot3D_TreeAndFields.m cannot be executed with current EstablishedLinks.dat file. ***\n');
    return
else
    fprintf('\n*** Executing Plot3D_TreeAndFields.m script. ***\n');
end

% Derive main parameters (number of nodes, spacings, and domain size):
N.x = Nxyz(1);        N.y = Nxyz(2);        N.z = Nxyz(3);
d.x = dxyz(1);        d.y = dxyz(2);        d.z = dxyz(3);      % in meters
L.x = (N.x-1)*d.x;    L.y = (N.y-1)*d.y;    L.z = (N.z-1)*d.z;  % in meters

% Initialize domain spaces for plotting:
fieldFactor = checkMagnitude(MaximumEfield);
spatialFactor = checkMagnitude([((0:(N.x-1))*d.x) ((0:(N.y-1))*d.y) ((0:(N.z-1))*d.z)]);
x = ((0:(N.x-1))*d.x)*spatialFactor.Number;
y = ((0:(N.y-1))*d.y)*spatialFactor.Number;
z = ((0:(N.z-1))*d.z + gnd.alt)*spatialFactor.Number;
[X,Y,Z] = meshgrid(x,y,z);
if strcmp(is.horizontal,'N')
    [vals.y,vals.z1] = meshgrid(y,z);
    [vals.x,vals.z2] = meshgrid(x,z);
elseif strcmp(is.horizontal,'Y')
    [vals.x,vals.y]  = meshgrid(x,y);
end

is.finalStep = size(Links.ID,1);

% Assigns plot height and width based on domain:
if ~isfield(sims,'plotWidth') || ~isfield(sims,'plotHeight')
    sims.plotWidth = 800;
    sims.plotHeight = round((Nxyz(3)/max(Nxyz(1:2)))*7)*80;
end

% Map ColorScale
gnd.color = [.75 .75 .75];
color     = zeros(is.finalStep,3);
linestyle = repelem("-",is.finalStep);
if strcmp(is.monoChrome,'N')
    for ii=1:is.finalStep
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

% File naming conventions:
if strcmp(is.whichScalar,'C')
    scalarType = 'Charge';
elseif strcmp(is.whichScalar,'P')
    scalarType = 'Potential';
end
if strcmp(is.zoomedIn,'Y')
    viewType = 'Zoomed';
elseif strcmp(is.zoomedIn,'N')
    viewType = 'FullScale';
end
if strcmp(is.horizontal,'Y')
    planeType = strcat('Horizontal-',num2str(floor(min(is.altitudes)*spatialFactor.Number)),spatialFactor.Unit,'m-',num2str(floor(max(is.altitudes)*spatialFactor.Number)),spatialFactor.Unit,'m');
elseif strcmp(is.horizontal,'N')
    planeType = 'Vertical';
end

% Set movie recording information:
if strcmp(is.Rec,'Y')
    Movie = VideoWriter(strcat(sims.pathVideos,'/Fields_',scalarType,'_',planeType,'_',viewType,'_',sims.objectName,'_',sims.objectType),'MPEG-4');
    if is.finalStep <= 60
        Movie.FrameRate = 4;
    elseif is.finalStep > 60 && is.finalStep <= 1000
        Movie.FrameRate = round(is.finalStep/20);
    else
        Movie.FrameRate = round(is.finalStep/60);
    end
    Movie.Quality = 100;
    open(Movie);
end

% Prepares the tiled layout for the figure:
if strcmp(is.highRes,'Y')
    is.resFactor = 2;
else
    is.resFactor = 1;
end
if strcmp(is.horizontal,'Y')
    figure(1); clf;
    tiledlayout(2,3,"TileSpacing","tight","Padding","tight")
    nexttile
    box on
    nexttile([2 2])
    hold on;
    set(gcf,'Position',[0,0,is.resFactor*4*round((1.8*sims.plotWidth)/3),is.resFactor*4*round((1.3*sims.plotHeight)/3)],'Color','White');
elseif strcmp(is.horizontal,'N')
    figure(2); clf;
    tiledlayout(1,3,"TileSpacing","tight","Padding","tight")
    nexttile(2)
    hold on;
    set(gcf,'Position',[0,0,is.resFactor*4*round((2.1*sims.plotWidth)/3),is.resFactor*4*round((sims.plotHeight)/3)],'Color','White');
end
set(gcf,'Resize','off');
grid on;
axis equal
if strcmp(is.zoomedIn,'Y')
    xlim(spatialFactor.Number*is.zoomedX);
    ylim(spatialFactor.Number*is.zoomedY);
    zlim(spatialFactor.Number*is.zoomedZ);
    arrowDensity = 2*round(10*(max(z)/(spatialFactor.Number*(is.zoomedZ(2)-is.zoomedZ(1)))));
    initiationMarker = 20*round(max(z)/(spatialFactor.Number*(is.zoomedZ(2)-is.zoomedZ(1))));
else
    xlim([0 max(x)]);
    ylim([0 max(y)]);
    zlim([0 max(z)]);
    arrowDensity = 10;
    if strcmp(is.horizontal,'N')
        initiationMarker = 20*is.resFactor;
        arrowDensity = 10*is.resFactor;
    elseif strcmp(is.horizontal,'Y')
        initiationMarker = 30*is.resFactor;
        arrowDensity = 7*is.resFactor;
    end
end

% Represents the neutrally charged (grounded) surface, if relevant:
if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
    P.x = [L.x 0 0 L.x]*spatialFactor.Number;
    P.y = [L.y L.y 0 0]*spatialFactor.Number;
    P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*spatialFactor.Number;
    patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
    if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
        patch(P.x, P.y, [z(end) z(end) z(end) z(end)], [z(end) z(end) z(end) z(end)],'FaceColor',gnd.color,'HandleVisibility','off');
        if strcmp(sims.BCtype,'TIN_CAN')
            patch([L.x L.x L.x L.x]*spatialFactor.Number, [L.y L.y 0 0]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
            patch([L.x L.x 0 0]*spatialFactor.Number, [L.y L.y L.y L.y]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
        end
    end
end
hold off
nexttile  
box on

% Initialize distance traveled for lightning links:
distance = 0;

% Converts altitudes to respective indices:  
is.altitudes = ((is.altitudes./d.z)+1); 

% For-loop to plot lightning discharge links and field lines:
for ii=0:is.finalStep
    if ii == 0 % (ambient case)
        nexttile(2)
        hold on
        if strcmp(is.whichScalar,'C')
            customTransparency = 0.25;
            scalar.data = load('../results/rhoAmb.dat','-ascii');
        elseif strcmp(is.whichScalar,'P')
            customTransparency = 0.04;
            scalar.data = load('../results/phiAmb.dat','-ascii');
        end
        scalar.data = ConvertTo3d(scalar.data,Nxyz); % nC/_m^3 or V
        scalar.max = max(max(max(scalar.data)));
        scalar.min = min(min(min(scalar.data)));
        if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
            patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
            if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                patch(P.x, P.y, [z(end) z(end) z(end) z(end)], [z(end) z(end) z(end) z(end)],'FaceColor',gnd.color,'HandleVisibility','off');
                if strcmp(sims.BCtype,'TIN_CAN')
                    patch([L.x L.x L.x L.x]*spatialFactor.Number, [L.y L.y 0 0]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
                    patch([L.x L.x 0 0]*spatialFactor.Number, [L.y L.y L.y L.y]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
                end 
            end
        end
        scalarValueUpdate('scalar',customTransparency,scalar,X,Y,Z,is,spatialFactor);
        x1 = Links.ID(1,1)*d.x*spatialFactor.Number;
        y1 = Links.ID(1,2)*d.y*spatialFactor.Number;
        z1 = ((Links.ID(1,3)*d.z)*spatialFactor.Number)+(gnd.alt*spatialFactor.Number);
        scatter3(x1,y1,z1,initiationMarker,'k','filled');
        titleformat = title(strcat(sims.disType," ",sims.objectType," on ",sims.objectName,": Ambient Conditions"),'FontSize',is.resFactor*20,'FontWeight','bold','Interpreter','latex');
        hold off

        % Determines the 2D electric field values for the step:
        E = consolidateEfield((ii-1),is,Nxyz);
        
        % Plots the electric field lines and ensures proper formatting for the first planar extract:
        nexttile(1)
        hold on
        if strcmp(is.horizontal,'Y')
            EfieldMag = sqrt(((E.x2Dmax).^2)+((E.y2Dmax).^2)+((E.z2Dmax).^2));
            direction = 2*atan(E.z2Dmax./sqrt(((E.x2Dmax).^2)+((E.y2Dmax).^2)))./pi;
        elseif strcmp(is.horizontal,'N')
            EfieldMag = sqrt(((E.y2D).^2)+((E.y2Dx).^2)+((E.z2Dy).^2));
            direction = 2*atan(E.z2Dy./sqrt(((E.y2Dx).^2)+((E.y2D).^2)))./pi;
        end
        direction(isnan(direction))=0;
        newMag = EfieldMag./MaximumEfield;
        tester = colormap(nexttile(1),createCustomColorMap('field',1));
        customColorData = zeros([size(direction,1),size(direction,2),3]);
        for i=1:size(direction,1)
            for j = 1:size(direction,2)
                for k = 1:3
                    customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                end
            end
        end
        p1.bar = colorbar;
        p1.bar.Ticks = [0, 0.25, 0.5, 0.75, 1];
        p1.bar.TickLabels = {'$E_z = -|\vec{E}|$','$E_z = -E_\parallel$','$E_z = 0$','$E_z = E_\parallel$','$E_z = |\vec{E}|$'};
        p1.bar.Label.Interpreter = 'latex';
        p1.bar.Label.String = strcat("Ratio of $E_z$ to $E_\parallel$ where $|\vec{E}|_\mathrm{max}$ = ",num2str(MaximumEfield*fieldFactor.Number,"%.1f")," ",fieldFactor.LaTeX,"V/m");
        p1.bar.Label.HorizontalAlignment = 'center';
        p1.bar.TickDirection = 'out';
        p1.bar.TickLabelInterpreter = 'latex';
        p1.bar.AxisLocation = "out";
        set(nexttile(1),'CLim',[0 1])
        set(nexttile(1),'CLimMode','manual')
        if strcmp(is.horizontal,'N')
            p1.bar.Label.FontSize = is.resFactor*14;
            p1.bar.Location = "southoutside"; 
            p1.mag = surf(vals.y,vals.z1,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
            p1.dir = streamslice(vals.y,vals.z1,E.y2D,E.z2Dy,arrowDensity,'arrows');
            set(gca,'FontSize',is.resFactor*10,'TickLabelInterpreter','latex');
            xlabel(strcat('$y$ (',spatialFactor.Unit,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            ylabel(strcat('$z$ (',spatialFactor.Unit,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            title(strcat("Vertical Electric Fields at $x$ = ",num2str(round((L.y*spatialFactor.Number)/2))," ",spatialFactor.LaTeX,"m"),'Interpreter','latex','FontSize',is.resFactor*20);
            axis xy
            axis image
            if strcmp(is.zoomedIn,'Y') 
                set(gca,'XDir','reverse','XLim',spatialFactor.Number*is.zoomedY,'YLim',spatialFactor.Number*is.zoomedZ);
            else
                set(gca,'XDir','reverse','XLim',[min(vals.y(:)) max(vals.y(:))],'YLim',[min(vals.z1(:)) max(vals.z1(:))]);
            end
        elseif strcmp(is.horizontal,'Y')
            p1.bar.Label.FontSize = is.resFactor*14;
            p1.bar.Location = "eastoutside"; 
            p1.mag = surf(vals.x,vals.y,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
            p1.dir = streamslice(vals.x,vals.y,E.x2Dmax,E.y2Dmax,arrowDensity,'arrows');
            set(gca,'FontSize',is.resFactor*10,'TickLabelInterpreter','latex');
            xlabel(strcat('$x$ (',spatialFactor.Unit,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            ylabel(strcat('$y$ (',spatialFactor.Unit,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            title(strcat("Horizontal Electric Fields at $z$ = ",num2str(z(max(is.altitudes)))," ",spatialFactor.LaTeX,"m"),'Interpreter','latex','FontSize',is.resFactor*20);
            axis xy
            axis image
            if strcmp(is.zoomedIn,'Y')
                set(gca,'XLim',spatialFactor.Number*is.zoomedX,'YLim',spatialFactor.Number*is.zoomedY);
            else
                set(gca,'XLim',[min(vals.x(:)) max(vals.x(:))],'YLim',[min(vals.y(:)) max(vals.y(:))]);
            end
        end
        box on
        hold off
       
        % Plots the electric field lines and ensures proper formatting for the second planar extract:
        if strcmp(is.horizontal,'N')
            nexttile(3)
            hold on
            EfieldMag = sqrt(((E.x2D).^2)+((E.x2Dy).^2)+((E.z2Dx).^2));
            direction = 2*atan(E.z2Dx./sqrt(((E.x2D).^2)+((E.x2Dy).^2)))./pi;
            tester = colormap(nexttile(3),createCustomColorMap('field',1));
        elseif strcmp(is.horizontal,'Y')
            nexttile(4)
            hold on
            EfieldMag = sqrt(((E.x2Dmin).^2)+((E.y2Dmin).^2)+((E.z2Dmin).^2));
            direction = 2*atan(E.z2Dmin./sqrt(((E.x2Dmin).^2)+((E.y2Dmin).^2)))./pi;
            tester = colormap(nexttile(4),createCustomColorMap('field',1));
        end
        direction(isnan(direction))=0;
        newMag = EfieldMag./MaximumEfield;
        customColorData = zeros([size(direction,1),size(direction,2),3]);
        for i=1:size(direction,1)
            for j = 1:size(direction,2)
                for k = 1:3
                    customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                end
            end
        end
        p2.bar = colorbar;
        p2.bar.Ticks = [0, 0.25, 0.5, 0.75, 1];
        p2.bar.TickLabels = {'$E_z = -|\vec{E}|$','$E_z = -E_\parallel$','$E_z = 0$','$E_z = E_\parallel$','$E_z = |\vec{E}|$'};
        p2.bar.Label.Interpreter = 'latex';
        p2.bar.Label.String = strcat("Ratio of $E_z$ to $E_\parallel$ where $|\vec{E}|_\mathrm{max}$ = ",num2str(MaximumEfield*fieldFactor.Number,"%.1f")," ",fieldFactor.LaTeX,"V/m");
        p2.bar.Label.HorizontalAlignment = 'center';
        p2.bar.TickDirection = 'out';
        p2.bar.TickLabelInterpreter = 'latex';
        p2.bar.AxisLocation = "out";
        if strcmp(is.horizontal,'N')
            p2.bar.Label.FontSize = is.resFactor*14;
            set(nexttile(3),'CLim',[0 1])
            set(nexttile(3),'CLimMode','manual')
            p2.bar.Location = "southoutside"; 
            p2.mag = surf(vals.x,vals.z2,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
            p2.dir = streamslice(vals.x,vals.z2,E.x2D,E.z2Dx,arrowDensity,'arrows');
            set(findobj('Type','line'),'Color',[0 0 0 0.5],'LineWidth',0.25);
            set(gca,'FontSize',is.resFactor*10,'TickLabelInterpreter','latex');
            xlabel(strcat('$x$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            ylabel(strcat('$z$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            title(strcat("Vertical Electric Fields at $y$ = ",num2str(round((L.x*spatialFactor.Number)/2))," ",spatialFactor.LaTeX,"m"),'Interpreter','latex','FontSize',is.resFactor*20);
            axis xy
            axis image
            if strcmp(is.zoomedIn,'Y')
                set(gca,'XLim',spatialFactor.Number*is.zoomedX,'YLim',spatialFactor.Number*is.zoomedZ);
            else
                set(gca,'XLim',[min(vals.x(:)) max(vals.x(:))],'YLim',[min(vals.z2(:)) max(vals.z2(:))]);
            end
        elseif strcmp(is.horizontal,'Y')
            p2.bar.Label.FontSize = is.resFactor*14;
            set(nexttile(4),'CLim',[0 1])
            set(nexttile(4),'CLimMode','manual')
            p2.bar.Location = "eastoutside"; 
            p2.mag = surf(vals.x,vals.y,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
            p2.dir = streamslice(vals.x,vals.y,E.x2Dmin,E.y2Dmin,arrowDensity,'arrows');
            set(findobj('Type','line'),'Color',[0 0 0 0.5],'LineWidth',0.25);
            set(gca,'FontSize',is.resFactor*10,'TickLabelInterpreter','latex');
            xlabel(strcat('$x$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            ylabel(strcat('$y$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',is.resFactor*16);
            title(strcat("Horizontal Electric Fields at $z$ = ",num2str(z(min(is.altitudes)))," ",spatialFactor.LaTeX,"m"),'Interpreter','latex','FontSize',is.resFactor*20);
            axis xy
            axis image
            if strcmp(is.zoomedIn,'Y')
                set(gca,'XLim',spatialFactor.Number*is.zoomedX,'YLim',spatialFactor.Number*is.zoomedY);                 
            else
                set(gca,'XLim',[min(vals.x(:)) max(vals.x(:))],'YLim',[min(vals.y(:)) max(vals.y(:))]);                                          
            end
            
        end
        clear EfieldMag; clear direction; clear newMag;
        box on
        hold off
        exportgraphics(gcf,strcat(sims.pathPNGs,'/Fields_',scalarType,'_',planeType,'_',viewType,'_Ambient_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',300);
    else % (all step values associated with the simulation)
        % Determines the 2D electric field values for the step:
        if strcmp(is.updateScalar,'Y') && (mod((ii-1),stepsaves) == 0 && (ismember(ii,is.setExtracts) || strcmp(is.Rec,'Y') || ismember(ii-1,stepsaves*floor((is.setExtracts-1)/stepsaves)))) || ii == is.finalStep
            E = consolidateEfield((ii-1),is,Nxyz);
            if strcmp(is.horizontal,'N')
                % Updates the electric field lines for the first planar extract:
                nexttile(1)
                delete(p1.mag); delete(p1.dir);
                hold on
                EfieldMag = sqrt(((E.y2D).^2)+((E.y2Dx).^2)+((E.z2Dy).^2));
                direction = 2*atan(E.z2Dy./sqrt(((E.y2Dx).^2)+((E.y2D).^2)))./pi;
                direction(isnan(direction))=0;
                newMag = EfieldMag./MaximumEfield;
                customColorData = zeros([size(direction,1),size(direction,2),3]);
                for i=1:size(direction,1)
                    for j = 1:size(direction,2)
                        for k = 1:3
                            customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                        end
                    end
                end
                p1.mag = surf(vals.y,vals.z1,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
                p1.dir = streamslice(vals.y,vals.z1,E.y2D,E.z2Dy,arrowDensity,'arrows');
                
                % Updates the electric field lines for the second planar extract:
                nexttile(3)
                delete(p2.mag); delete(p2.dir);
                hold on
                EfieldMag = sqrt(((E.x2D).^2)+((E.x2Dy).^2)+((E.z2Dx).^2));
                direction = 2*atan(E.z2Dx./sqrt(((E.x2D).^2)+((E.x2Dy).^2)))./pi;
                direction(isnan(direction))=0;
                newMag = EfieldMag./MaximumEfield;
                customColorData = zeros([size(direction,1),size(direction,2),3]);
                for i=1:size(direction,1)
                    for j = 1:size(direction,2)
                        for k = 1:3
                            customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                        end
                    end
                end
                p2.mag = surf(vals.x,vals.z2,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
                p2.dir = streamslice(vals.x,vals.z2,E.x2D,E.z2Dx,arrowDensity,'arrows');
                hold off
            elseif strcmp(is.horizontal,'Y')
                % Updates the electric field lines for the first planar extract:
                nexttile(1)
                delete(p1.mag); delete(p1.dir);
                hold on
                EfieldMag = sqrt(((E.x2Dmax).^2)+((E.y2Dmax).^2)+((E.z2Dmax).^2));
                direction = 2*atan(E.z2Dmax./sqrt(((E.x2Dmax).^2)+((E.y2Dmax).^2)))./pi;
                direction(isnan(direction))=0;
                newMag = EfieldMag./MaximumEfield;
                customColorData = zeros([size(direction,1),size(direction,2),3]);
                for i=1:size(direction,1)
                    for j = 1:size(direction,2)
                        for k = 1:3
                            customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                        end
                    end
                end
                p1.mag = surf(vals.x,vals.y,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
                p1.dir = streamslice(vals.x,vals.y,E.x2Dmax,E.y2Dmax,arrowDensity,'arrows');
                
                % Updates the electric field lines for the second planar extract:
                nexttile(4)
                delete(p2.mag); delete(p2.dir);
                hold on
                EfieldMag = sqrt(((E.x2Dmin).^2)+((E.y2Dmin).^2)+((E.z2Dmin).^2));
                direction = 2*atan(E.z2Dmin./sqrt(((E.x2Dmin).^2)+((E.y2Dmin).^2)))./pi;
                direction(isnan(direction))=0;
                newMag = EfieldMag./MaximumEfield;
                customColorData = zeros([size(direction,1),size(direction,2),3]);
                for i=1:size(direction,1)
                    for j = 1:size(direction,2)
                        for k = 1:3
                            customColorData(i,j,k) = ((1-newMag(i,j))*(1-tester(round((1+direction(i,j))*50)+1,k)))+tester(round((1+direction(i,j))*50)+1,k);
                        end
                    end
                end
                p2.mag = surf(vals.x,vals.y,-1+zeros(size(vals.x)),customColorData,'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
                p2.dir = streamslice(vals.x,vals.y,E.x2Dmin,E.y2Dmin,arrowDensity,'arrows');
                hold off
            end
            set(findobj('Type','line'),'Color',[0 0 0 0.5],'LineWidth',0.25);
            clear EfieldMag; clear direction; clear newMag;
        
            % Updates charge density/potential plot:
            nexttile(2)
            hold on
            allPatches = findall(nexttile(2),'type','patch');
            delete(allPatches);
            if strcmp(is.whichScalar,'C')
                if ii == is.finalStep
                    scalar.data = load('../results/rho3d.dat','-ascii');
                else
                    scalar.data = load(['../results/rho3d',num2str(ii-1),'.dat'],'-ascii');
                end
            elseif strcmp(is.whichScalar,'P')
                if ii == is.finalStep
                    scalar.data = load('../results/phi.dat','-ascii');
                else
                    scalar.data = load(['../results/phi3d',num2str(ii-1),'.dat'],'-ascii');
                end
            end
        
            scalar.data = ConvertTo3d(scalar.data,Nxyz); % nC/_m^3 or V
            if scalar.max < max(max(max(scalar.data))) || scalar.max > 10*max(max(max(scalar.data)))
                scalar.max = max(max(max(scalar.data)));
            end
            if abs(scalar.min) < abs(min(min(min(scalar.data)))) || abs(scalar.min) > 10*abs(min(min(min(scalar.data))))
                scalar.min = min(min(min(scalar.data)));
            end
            if strcmp(sims.BCtype,'G') || strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
                if strcmp(sims.BCtype,'G_G') || strcmp(sims.BCtype,'TIN_CAN')
                    patch(P.x, P.y, [z(end) z(end) z(end) z(end)], [z(end) z(end) z(end) z(end)],'FaceColor',gnd.color,'HandleVisibility','off'); 
                    if strcmp(sims.BCtype,'TIN_CAN')
                        patch([L.x L.x L.x L.x]*spatialFactor.Number, [L.y L.y 0 0]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
                        patch([L.x L.x 0 0]*spatialFactor.Number, [L.y L.y L.y L.y]*spatialFactor.Number, [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number], [gnd.alt*spatialFactor.Number z(end) z(end) gnd.alt*spatialFactor.Number],'FaceColor',gnd.color,'HandleVisibility','off'); 
                    end
                end
            end
            scalarValueUpdate('scalar',customTransparency,scalar,X,Y,Z,is,spatialFactor);
            hold off
        end    
          
        % Revisiting panel 2 to draw the link propagations:
        nexttile(2)
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
        if strcmp(is.horizontal,'N')
            plot3([x1, x2]*spatialFactor.Number,[y1, y2]*spatialFactor.Number,[z1, z2]*spatialFactor.Number,'Color',color(ii,:),'LineStyle',linestyle(ii),'HandleVisibility','off','LineWidth',1.5);
        elseif strcmp(is.horizontal,'Y')
            plot3([x1, x2]*spatialFactor.Number,[y1, y2]*spatialFactor.Number,[z1, z2]*spatialFactor.Number,'Color',color(ii,:),'LineStyle',linestyle(ii),'HandleVisibility','off','LineWidth',2);
        end
    end
    % Formatting title to reflect step value:
    nexttile(2)
    delete(titleformat);
    titleformat = title(strcat(sims.disType," ",sims.objectType," on ",sims.objectName," after ", int2str(ii) ," step(s)"),'FontSize',is.resFactor*20,'FontWeight','bold','Interpreter','latex');
    box off
    if strcmp(is.singleImages,'Y') || ismember(ii,is.setExtracts)
        exportgraphics(gcf,strcat(sims.pathPNGs,'/Fields_',scalarType,'_',planeType,'_',viewType,'_Step-',num2str(ii),'_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',300);
    end
    if strcmp(is.Rec,'Y')
        frame = getframe(gcf);
        writeVideo(Movie,frame);
    end
end

% Outputs general information about the discharge's propagation:
fprintf(strcat('\n\t',sims.objectType," reaches minimum of ",num2str(minHeight*spatialFactor.Number,'%.2f')," ",spatialFactor.Prefix,'meters (',num2str((minHeight-initHeight)*spatialFactor.Number,'%+.2f')," ",spatialFactor.Prefix,'meters)\n'));
fprintf(strcat('\n\t',sims.objectType," initiated at ",num2str(initHeight*spatialFactor.Number,'%.2f')," ",spatialFactor.Prefix,'meters\n'));
fprintf(strcat('\n\t',sims.objectType," reaches maximum of ",num2str(maxHeight*spatialFactor.Number,'%.2f')," ",spatialFactor.Prefix,'meters (',num2str((maxHeight-initHeight)*spatialFactor.Number,'%+.2f')," ",spatialFactor.Prefix,'meters)\n'));
fprintf(strcat('\n\t',sims.objectType," has propagated a total of ",num2str(distance*spatialFactor.Number,'%.2f')," ",spatialFactor.Prefix,'meters\n'));
hold off;

% Finalizes the movie, if relevant:
if strcmp(is.Rec,'Y')
    frame = getframe(gcf);
    writeVideo(Movie,frame);
    close(Movie);
end

% Exports the final graphic:
delete(titleformat);
titleformat = title(strcat(sims.disType," ",sims.objectType," on ",sims.objectName,": Final Results"),'FontSize',is.resFactor*20,'FontWeight','bold','Interpreter','latex');
exportgraphics(gcf,strcat(sims.pathPNGs,'/Fields_',scalarType,'_',planeType,'_',viewType,'_Final_',sims.objectName,'_',sims.objectType,'.png'),'Resolution',300);

%% Functions:
% Converts 3D electric field data files into respective 2D components:
function E = consolidateEfield(num,is,Nxyz)
    if num == -1 % (ambient case)
        Ex3D     = ConvertTo3d(load('../results/Ex3dAmb.dat'),Nxyz);
        Ey3D     = ConvertTo3d(load('../results/Ey3dAmb.dat'),Nxyz);
        Ez3D     = ConvertTo3d(load('../results/Ez3dAmb.dat'),Nxyz);
    elseif num+1 == is.finalStep % (final result)
        Ex3D     = ConvertTo3d(load('../results/Ex3d.dat'),Nxyz);
        Ey3D     = ConvertTo3d(load('../results/Ey3d.dat'),Nxyz);
        Ez3D     = ConvertTo3d(load('../results/Ez3d.dat'),Nxyz);
    else % (all other steps)
        Ex3D     = ConvertTo3d(load(['../results/Ex3d',num2str(num),'.dat']),Nxyz);
        Ey3D     = ConvertTo3d(load(['../results/Ey3d',num2str(num),'.dat']),Nxyz);
        Ez3D     = ConvertTo3d(load(['../results/Ez3d',num2str(num),'.dat']),Nxyz);
    end
    if strcmp(is.horizontal,'Y')
        E.x2Dmin = permute(Ex3D(:,:,min(is.altitudes)),[2,1,3]);
        E.x2Dmax = permute(Ex3D(:,:,max(is.altitudes)),[2,1,3]);
        E.y2Dmin = permute(Ey3D(:,:,min(is.altitudes)),[2,1,3]);
        E.y2Dmax = permute(Ey3D(:,:,max(is.altitudes)),[2,1,3]);
        E.z2Dmin = permute(Ez3D(:,:,min(is.altitudes)),[2,1,3]);
        E.z2Dmax = permute(Ez3D(:,:,max(is.altitudes)),[2,1,3]);
    elseif strcmp(is.horizontal,'N')
        E.x2D    = permute(Ex3D(:,round((Nxyz(2)-1)/2),:),[3,1,2]);
        E.x2Dy   = permute(Ey3D(:,round((Nxyz(2)-1)/2),:),[3,1,2]);
        E.y2D    = permute(Ey3D(round((Nxyz(1)-1)/2),:,:),[3,2,1]);
        E.y2Dx   = permute(Ex3D(round((Nxyz(1)-1)/2),:,:),[3,2,1]);
        E.z2Dx   = permute(Ez3D(:,round((Nxyz(2)-1)/2),:),[3,1,2]);
        E.z2Dy   = permute(Ez3D(round((Nxyz(1)-1)/2),:,:),[3,2,1]);
    end
end

% Similar functionality to plottingChargeRegions.m file, but expanded for 
% all scalar quantities plus a custom colorbar placement for a compact form:
function scalarValueUpdate(colorbarRange,alphaValue,scalarvaluesOG,Xval,Yval,Zval,is,spatialFactor)
    % Determines the range of the colorbar among other factors:
    if strcmp(is.whichScalar,'C') 
        scalarFactor = checkMagnitude((10^(-9)*scalarvaluesOG.data(:)));
        scalarFactor.Number = (scalarFactor.Number)*(10^(-9));
    elseif strcmp(is.whichScalar,'P')
        scalarFactor = checkMagnitude(scalarvaluesOG.data(:));
    end
    tol = log10(round(max(max(max(abs(scalarvaluesOG.data)))),2,'significant'));   
    if tol < 0 
        tol = floor(tol);
    else
        tol = ceil(tol);
    end
    scalarvalues.data = (scalarFactor.Number)*round(scalarvaluesOG.data,(2-tol));
    uniqueValues = unique(nonzeros(scalarvalues.data));
    
    % Ensures values are plotted from weakest to strongest regardless of polarity:
    testFlag = 1;
    while testFlag < length(uniqueValues) && uniqueValues(testFlag+1)<0
        testFlag = testFlag+1;
    end
    uniqueValues(1:testFlag) = uniqueValues(testFlag:-1:1);
    scalarvalues.max = max([scalarFactor.Number*scalarvaluesOG.max max(uniqueValues)]);
    scalarvalues.min = min([scalarFactor.Number*scalarvaluesOG.min min(uniqueValues)]);

    % Prevents misread of charge density/potential when cloud height is less than dz:
    testingvalues = zeros([length(uniqueValues),1]);
    testinglocation = zeros([length(uniqueValues),1]);
    for testloop = 1:1:length(scalarvalues.data(1,1,:))
        testingunique = unique(nonzeros(scalarvalues.data(:,:,testloop)));
        if testingunique~=0
            [~, testloc] = ismember(testingunique,uniqueValues);
            testingvalues(testloc) = testingvalues(testloc)+1;
            testinglocation(testloc) = testloop;
        end
    end

    % Creates a unique colormap to represent the charge regions:
    rgbValues = createCustomColorMap(colorbarRange,1);
    rgbValuesAdjusted = createCustomColorMap(colorbarRange,alphaValue);
   
    % Copies data to a nearby height without overwriting pre-existing data:
    for testloop2 = 1:1:length(uniqueValues)
        if testingvalues(testloop2)==1
            if max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)-1))))==0 && max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)+1))))==0
                scalarvalues.data(:,:,testinglocation(testloop2)-1)=scalarvalues.data(:,:,testinglocation(testloop2));
            elseif max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)-1))))~=0 && max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)+1))))==0
                scalarvalues.data(:,:,testinglocation(testloop2)+1)=scalarvalues.data(:,:,testinglocation(testloop2));
            elseif max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)-1))))==0 && max(abs(unique(scalarvalues.data(:,:,testinglocation(testloop2)+1))))~=0
                scalarvalues.data(:,:,testinglocation(testloop2)-1)=scalarvalues.data(:,:,testinglocation(testloop2));
            end
        end
    end

    % Determines truly unique values for legend usage:
    trueUniqueValues = uniqueValues;
    middle = (((length(rgbValues(:,1))-1)/2)+1);
    midRange = floor((0.025*(length(rgbValues(:,1))-1)))-1;
    
    % Removes values that are approximately zero (i.e. neutral):
    for i = length(trueUniqueValues):-1:1
        [nullInd,~] = colorDetermination(trueUniqueValues(i),max(abs(uniqueValues)),rgbValues);
        if nullInd <= (middle+round(midRange/2)) && nullInd >= (middle-round(midRange/2))
            trueUniqueValues(i)=[];
        end
    end
    max_scalar_value = max(trueUniqueValues);
    min_scalar_value = min(trueUniqueValues);
    %fprintf(['Maximum value is ',num2str(max(trueUniqueValues)),'\nMinimum value is ',num2str(min(trueUniqueValues)),'\nLength of unique values is ',num2str(length(trueUniqueValues))]);
    
    % Plots isometric regions in a representative color:
    colorIndices = zeros([length(uniqueValues),1]);
    colorVertices = zeros([length(uniqueValues),3]);
    for j = length(uniqueValues):-1:1
        [colorIndices(j),colorVertices(j,:)] = colorDetermination(uniqueValues(j),max(abs(uniqueValues)),rgbValues);
        
        % If the region is nonzero:
        if colorIndices(j) ~= (((length(rgbValues(:,1))-1)/2)+1)
            % Is the current value a 'truly unique' value?
            [included, location] = ismember(uniqueValues(j),trueUniqueValues);
            % If so, add it to the legend:
            if included == 1
                if trueUniqueValues(location)>0
                    p0 = patch(isosurface(Xval,Yval,Zval,-1*pagetranspose(scalarvalues.data),-trueUniqueValues(location)));
                else
                    p0 = patch(isosurface(Xval,Yval,Zval,pagetranspose(scalarvalues.data),trueUniqueValues(location)));
                end
                if length(trueUniqueValues)>3
                    if colorIndices(j) <= (middle+midRange) && colorIndices(j) >= (middle-midRange)
                        chosenAlpha = alphaValue/3;
                    else
                        chosenAlpha = alphaValue;
                    end
                    if uniqueValues(location) == max_scalar_value && min_scalar_value > 0
                        set(p0,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','on','DisplayName',strcat('$$0 \leq \rho \leq$$ ',num2str(trueUniqueValues(location)),' ',scalarFactor.LaTeX,'C/m$^3$'));
                        isonormals(Xval,Yval,Zval,pagetranspose(scalarvalues.data),p0);
                        drawnow
                    elseif uniqueValues(location) == min_scalar_value && max_scalar_value < 0
                        set(p0,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','on','DisplayName',strcat('$$0 \geq \rho \geq$$ ',num2str(trueUniqueValues(location)),' ',scalarFactor.LaTeX,'C/m$^3$'));
                        isonormals(Xval,Yval,Zval,pagetranspose(scalarvalues.data),p0);
                        drawnow
                    else
                        set(p0,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',chosenAlpha,'HandleVisibility','off');
                        isonormals(Xval,Yval,Zval,pagetranspose(scalarvalues.data),p0);
                        drawnow
                    end
                else
                    set(p0,'FaceColor',colorVertices(j,:),'EdgeAlpha',0,'FaceAlpha',alphaValue,'HandleVisibility','on','DisplayName',strcat('Charge Density $$\approx$$ ',num2str(trueUniqueValues(location)),' ',scalarFactor.LaTeX,'C/m$^3$'));
                    isonormals(Xval,Yval,Zval,pagetranspose(scalarvalues.data),p0);
                    drawnow
                end
            end
        end
    end
    
    % Resolves formatting issues with the custom colorbar:
    cmap = (rgbValues+rgbValuesAdjusted+rgbValuesAdjusted)/3;
    clim([-max(abs([scalarvalues.min scalarvalues.max max(abs(uniqueValues))])) max(abs([scalarvalues.min scalarvalues.max max(abs(uniqueValues))]))]);
    colormap(gca,cmap);
    if size(findall(gcf,'type','colorbar'),1) == 0
        c = colorbar;
        c.Label.Interpreter = 'latex';
        if strcmp(is.whichScalar,'C')
            c.Label.String = strcat('Total Charge Density (',scalarFactor.LaTeX,'C/m$^3$)');
        elseif strcmp(is.whichScalar,'P')
            c.Label.String = strcat('Total Potential (',scalarFactor.LaTeX,'V)');
        end
        c.Label.HorizontalAlignment = 'center';
        c.Label.FontSize = 16*is.resFactor;
        c.TickDirection = 'out';
        c.TickLabelInterpreter = 'latex';
        c.AxisLocation = "out";
        if strcmp(is.horizontal,'N')
            c.Location = "southoutside"; 
        elseif strcmp(is.horizontal,'Y')
            c.Location = "eastoutside";
        end
    end
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 10*is.resFactor;
    xlabel(strcat('$x$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',16*is.resFactor,'HorizontalAlignment','left');
    ylabel(strcat('$y$ (',spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',16*is.resFactor,'HorizontalAlignment','right');
    grid on
    view(-45,5)
end