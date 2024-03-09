% Initiate
close all
clearvars -except sims is
clf

% For formatting exported image size:
%% Streamer
%positionWidth = 600;
%positionHeight = 1200;
% Leader
positionWidth = 800;
positionHeight = 800;

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

if isempty(Links.ID)
    fprintf('\n*** LightningVisual.m cannot be executed with current EstablishedLinks.dat file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing LightningVisual.m script. ***\n');
    % User-Based Settings:
    prompt_Rec = '\nWould you like to record the lightning propagation as a movie? (Y / N)\n-->';
    is.Rec = input(prompt_Rec,'s');                    
    while ~strcmp(is.Rec,'Y') && ~strcmp(is.Rec,'N')
        fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
        is.Rec = input(prompt_Rec,'s');
    end
    prompt_updateChargeDensity = '\nWould you like to update the charge distribution coloring for every saved timestep? (Y / N)\n-->';
    is.updateChargeDensity = input(prompt_updateChargeDensity,'s');                    
    while ~strcmp(is.updateChargeDensity,'Y') && ~strcmp(is.updateChargeDensity,'N')
        fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
        is.updateChargeDensity = input(prompt_updateChargeDensity,'s');
    end
    prompt_highResolution = '\nWould you like to save the final image as a very high resolution image? (Y / N)\nWARNING: Only recommended for preparing posters.\n-->';
    is.highResolution = input(prompt_highResolution,'s');                    
    while ~strcmp(is.highResolution,'Y') && ~strcmp(is.highResolution,'N')
        fprintf('\n\tNot an acceptable input. Please enter Y (for yes) or N (for no).\n');
        is.highResolution = input(prompt_highResolution,'s');
    end
    prompt_Grounded = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
    is.Grounded = input(prompt_Grounded,'s');                    
    while ~strcmp(is.Grounded,'FS') && ~strcmp(is.Grounded,'G')
        fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
        is.Grounded = input(prompt_Grounded,'s');
    end
end

% Derive main parameters
Links.Nb = size(Links.ID);
Links.Nb = Links.Nb(1);
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

clear dxyz
clear InitPoint
cd ../viz

rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3

% Linear spaces for the three position dimensions:
x = ((0:(N.x-1))*d.x)*1e-3;
y = ((0:(N.y-1))*d.y)*1e-3;
z = ((0:(N.z-1))*d.z + gnd.alt)*1e-3;

[X,Y,Z]=meshgrid(x,y,z);
rho.max = .95* max(max(max(rho.data)));
rho.min = .95* min(min(min(rho.data)));

% Map ColorScale
gnd.color = [.75 .75 .75];%[.718 .255 .055];
color     = colormap(jet(Links.Nb));
is.BW = 1; %input('Is plot Monochrome? (1: yes, else: no)\n>> ');
if (is.BW == 1)
    for ii=1:Links.Nb
        color(ii,:) = [0 0 0];
    end
end
% Set movie recording
if (strcmp(is.Rec,'Y') == 1)
    Movie = VideoWriter([sims.pathVideos,'/',sims.objectName,'_',sims.objectType,'Video'],'MPEG-4');
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
figure(1);
set(gcf,'Position',[0,0,positionWidth,positionHeight]);
set(gcf,'Resize','off')
hold on;
grid on;

% Sets bounds for the axes (comment out if clouds get cut off):
axis([L.x*1/5 L.x*4/5 L.y*1/5 L.y*4/5 gnd.alt 2/2*(L.z+gnd.alt)]*1e-3) % Slight crop
%axis([0 L.x 0 L.y gnd.alt 2/2*(L.z+gnd.alt)]*1e-3)                     % Full span 

% Plots the cloud structure with the defined function below:
axis equal
xlim([0 max(x)]);
ylim([0 max(y)]);
%zlim([25 70]); % to crop the altitude range for visualization
zlim([0 max(z)]);
set(legend,'Position',[0.225 0.7 .5 .0375],'box','off')
%set(legend,'location','southoutside','box','on')
set(gcf,'Resize','off')

% Represents the neutrally charged (grounded) surface:
if strcmp(is.Grounded,'G')
    P.x = [L.x 0 0 L.x]*1e-3;
    P.y = [L.y L.y 0 0]*1e-3;
    P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
    patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');
end

% Initialize distance traveled for lightning links:
distance = 0;
% For-loop to plot lightning discharge links:
for ii=1:Links.Nb
    %     if(rem(ii,10)==0)
    %         %         pause;
    %     end
    if mod((ii-1),stepsaves) == 0
        if ii == 1
            plottingChargeRegions('white',0.25,rho,X,Y,Z);
            pause
        else
            if strcmp(is.updateChargeDensity,'Y') == 1
                allPatches = findall(gcf,'type','patch');
                delete(allPatches);
                rho.data = load(['../results/rho3d',num2str(ii-1),'.dat'],           '-ascii');
                rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3
                rho.max = .95* max(max(max(rho.data)));
                rho.min = .95* min(min(min(rho.data)));
                plottingChargeRegions('white',0.4,rho,X,Y,Z);
            end
        end
    end    
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
    
    % Formatting axes:
    set(gcf,'Position',[0,0,positionWidth,positionHeight]);
    set(gcf,'Resize','off')
    %axis equal
    %{
    xticks([0 5 10 15 20]);
    xticklabels({'0','5','10','12'});
    yticks([0 4 8 12]);
    yticklabels({'0','4','8','12'});
    zticks([46 50 54 58 62 66 70]);
    zticklabels({'46','50','54','58','62','66','70'});
    %}
    box on
    title([sims.objectType,' discharge after ', int2str(ii) ,' step(s)'],'FontSize',28,'FontWeight','bold','Interpreter','latex');
    if(strcmp(is.Rec,'Y') == 1)
        set(gcf,'Position',[0,0,positionWidth,positionHeight]); 
        set(gcf,'Resize','off')
        frame = getframe(gcf);
        writeVideo(Movie,frame);
    end
end
fprintf(['\n\t',sims.objectType,' reaches minimum of ',num2str(minHeight,'%.2f'),' meters (',num2str(minHeight-initHeight,'%+.2f'),' meters)\n']);
fprintf(['\n\t',sims.objectType,' initiated at ',num2str(initHeight,'%.2f'),' meters\n']);
fprintf(['\n\t',sims.objectType,' reaches maximum of ',num2str(maxHeight,'%.2f'),' meters (',num2str(maxHeight-initHeight,'%+.2f'),' meters)\n']);
fprintf(['\n\t',sims.objectType,' has propagated a total of ',num2str(distance,'%.2f'),' meters\n']);
%camlight; lighting gouraud

hold off;
% Record the movie
if (strcmp(is.Rec,'Y') == 1)
    frame = getframe(gcf);
    writeVideo(Movie,frame);
    close(Movie);
end
%title('(b)','Interpreter','latex','FontSize',32,'Units','normalized');
%titleInfo = get(gca,'title');
%set(titleInfo,'Position', [((titleInfo.Extent(3))-(((titleInfo.Parent.InnerPosition(1)/titleInfo.Parent.InnerPosition(3)))/titleInfo.Parent.Position(3))-(titleInfo.Parent.Position(1)/titleInfo.Parent.OuterPosition(3))) 1.004 0]);
title(['Simulated ', sims.objectType,' Discharge: ',sims.objectName],'FontSize',28,'FontWeight','bold','Interpreter','latex');    
set(gcf,'Position',[0,0,positionWidth,positionHeight]); 
set(gcf,'Resize','off')
% If the 'export_fig' function is assigned to the pathtool:
if exist('export_fig') == 2 && strcmp(is.highResolution,'Y') == 1
    currentFolder = pwd;
    cd(sims.pathPNGs);
    export_fig HighRes_Discharge.png -transparent -m8
    cd(currentFolder);
else
    exportgraphics(gcf,[sims.pathPNGs,'/Lightning_',sims.objectName,'_',sims.objectType,'.png'],'Resolution',600);
end