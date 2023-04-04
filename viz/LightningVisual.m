% Initiate
close all
clearvars
clf
clc

% User-Based Settings:
objectName = 'Jupiter';         % Name of planetary body being visualized
objectType = 'Streamer';        % Type of discharge simulated
is.Rec = 0;                     % Record the movie? (1: yes, else: no)
is.updateChargeDensity = 0;     % Update charge distribution coloring every saved timestep (1: yes, else: no)
is.highResolution = 0;          % will save image as 51M pixel image if set to 1, use for posters only!

% Settings to ensure proper directory referencing:
videoPath = ['../Figures/',objectName,'/',objectType,'/Videos'];
imagePath = ['../Figures/',objectName,'/',objectType,'/PNGs'];
if ~exist(videoPath,'dir')
    mkdir(videoPath);
end
if ~exist(imagePath,'dir')
    mkdir(imagePath);
end
cd ../results

% Load data files
load('dxyz.dat',             '-ascii');
load('Nxyz.dat',             '-ascii');
load('InitPoint.dat',        '-ascii');
rho.data = load('rhoAmb.dat',           '-ascii');
Links.ID = load('EstablishedLinks.dat', '-ascii');
gnd.alt  = load('z_gnd.dat',            '-ascii');

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

rho.data = ConvertTo3d(rho.data,Nxyz); % _C/_m^3

clear dxyz
clear InitPoint
cd ../viz

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
if (is.Rec == 1)
    Movie = VideoWriter([videoPath,'/',objectName,'_',objectType,'Video'],'MPEG-4');
    open(Movie);
end
% Draw the tree
figure(1);
set(gcf,'Position',[0,0,800,1000]);
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
P.x = [L.x 0 0 L.x]*1e-3;
P.y = [L.y L.y 0 0]*1e-3;
P.z = [gnd.alt gnd.alt gnd.alt gnd.alt]*1e-3;
patch(P.x, P.y, P.z, gnd.alt,'FaceColor',gnd.color,'HandleVisibility','off');

% Initialize distance traveled for lightning links:
distance = 0;
% For-loop to plot lightning discharge links:
for ii=1:Links.Nb
    %     if(rem(ii,10)==0)
    %         %         pause;
    %     end
    if mod((ii-1),100) == 0
        if ii == 1
            plottingChargeRegions('white',0.25,rho,X,Y,Z);
        else
            if is.updateChargeDensity == 1
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
    
    % Summing lightning link to overall distance of link traveled (m):
    distance = distance + sqrt(((x2-x1)^2) + ((y2-y1)^2) + ((z2-z1)^2));

    % Plotting link:
    plot3(...
        [x1, x2]*1e-3,...
        [y1, y2]*1e-3,...
        [z1, z2]*1e-3,...
        'Color',color(ii,:),'HandleVisibility','off');
    
    % Formatting axes:
    set(gcf,'Position',[0,0,800,1000]);
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
    title([objectType,' discharge after ', int2str(ii) ,' step(s)'],'FontSize',28,'FontWeight','bold','Interpreter','latex');
    if(is.Rec == 1)
        set(gcf,'Position',[0,0,800,1000]); 
        set(gcf,'Resize','off')
        frame = getframe(gcf);
        writeVideo(Movie,frame);
    end
end
fprintf(['\n',objectType,' has propagated %.2f meters\n',distance]);
pause
%camlight; lighting gouraud


hold off;
% Record the movie
if (is.Rec == 1)
    frame = getframe(gcf);
    writeVideo(Movie,frame);
    close(Movie);
end
title(['Simulated ',objectType,' Discharge: ',objectName],'FontSize',28,'FontWeight','bold','Interpreter','latex');
set(gcf,'Position',[0,0,800,1000]); 
set(gcf,'Resize','off')
% If the 'export_fig' function is assigned to the pathtool:
if exist('export_fig') == 2 && is.highResolution == 1
    highResSaveFile = [imagePath,'/',objectName,'_',objectType,'_HighRes.png'];
    export_fig ../Figures/HighRes_Discharge.png -transparent -m8
else
    exportgraphics(gcf,[imagePath,'/',objectName,'_',objectType,'.png'],'Resolution',600);
    exportgraphics(gcf,[imagePath,'/',objectName,'_',objectType,'.eps'],'Resolution',600);
end
    
function [AA] = ConvertTo3d(A,B)
    [M, N] = size(A);
    AA = zeros(B');
    for m=1:M
        for n=1:N
            ii = rem(m,B(1));
            if(ii==0)
                ii = B(1);
            end
            jj = n;
            kk = (m-ii)/B(1)+1;
            AA(ii,jj,kk) = A(m,n);
        end
    end
end