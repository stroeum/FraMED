clearvars -except sims

%% Check pathway linking:
if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
    prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
    sims.objectType = input(prompt2,'s');
    while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
        fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
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
if ~exist('summary.txt','file')
    fprintf('\n*** Plot1D_RuntimeResults.m cannot be executed with current summary.txt file. ***\n');
    cd ../viz
    return
else
    fprintf('\n*** Executing Plot1D_RuntimeResults.m script. ***\n');
end

%% Data analysis:
% Open summary textfile:
fid   = fopen('summary.txt','r');
tline = fgetl(fid);
cd ../viz

% Initialize number of links:
links.tracker       = 0;
while ischar(tline)
    tline = fgetl(fid);
    if size(tline,2)>27 && strcmp(tline(1:27),'ii:	 Number of candidates: ')
        links.tracker = links.tracker + 1;
        Candidates.values(links.tracker) = [str2double(tline(28:end))];
    elseif size(tline,2)>34 && strcmp(tline(1:34),'ii:	 Maximum candidate overreach: ')
        Overreach.values(links.tracker) = [str2double(tline(35:(end-1)))];
    elseif size(tline,2)>33 && strcmp(tline(1:33),'ii:	 Run time for Link addition: ')
        LinkRuntime.values(links.tracker) = [str2double(tline(34:(end-2)))];
    elseif size(tline,2)>33 && strcmp(tline(1:33),'ii:	 Run time for Qminimization: ')
        QMinimization.values(links.tracker) = [str2double(tline(34:(end-1)))];
    end
end
% Close file:
fclose(fid);

% Create step number array:
links.array = (1:1:links.tracker)';

% Determining min/max for Candidates:
[Candidates.max, Candidates.maxIndex] = max(Candidates.values);
[Candidates.min, Candidates.minIndex] = min(Candidates.values);

% Determining min/max for Overreach:
[Overreach.max, Overreach.maxIndex] = max(Overreach.values);
[Overreach.min, Overreach.minIndex] = min(Overreach.values);

% Determining min/max for LinkRuntime:
[LinkRuntime.max, LinkRuntime.maxIndex] = max(LinkRuntime.values);
[LinkRuntime.min, LinkRuntime.minIndex] = min(LinkRuntime.values);

% Determining min/max for QMinimization:
if exist('QMinimization','var')
    [QMinimization.max, QMinimization.maxIndex] = max(QMinimization.values);
    [QMinimization.min, QMinimization.minIndex] = min(QMinimization.values);
end

%% Plot Formatting:
% Initialize figure layout:
clf
if exist('QMinimization','var')
    tiledlayout(2,2, 'Padding', 'tight', 'TileSpacing', 'compact','Position',[0 0 800 600], 'PositionConstraint','outerposition'); 
else
    tiledlayout(1,3, 'Padding', 'tight', 'TileSpacing', 'loose','Position',[0 0 1200 400], 'PositionConstraint','outerposition');
end

% Plot number of candidates (row 1, column 1):
nexttile
hold on
plot(links.array,Candidates.values);
set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(Candidates.values))],'TickLabelInterpreter','latex');
xlabel('Step number','Interpreter','latex','FontSize',16);
ylabel('Number of candidate nodes','Interpreter','latex','FontSize',16);
legend(['Minimum candidates: ',num2str(Candidates.min),' nodes (step ',num2str(Candidates.minIndex),')',newline,'Maximum candidates: ',num2str(Candidates.max),' nodes (step ',num2str(Candidates.maxIndex),')'],'Interpreter','latex','location','best', 'box','off','FontSize',10)
box on
hold off

% Plot candidate overreach percentages (row 1, column 2):
nexttile
hold on
plot(links.array,Overreach.values);
if exist('QMinimization','var')
    set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(Overreach.values))],'TickLabelInterpreter','latex','YAxisLocation','right');
else
    set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(Overreach.values))],'TickLabelInterpreter','latex');
end
xlabel('Step number','Interpreter','latex','FontSize',16);
ylabel('Maximum candidate overreach (\%)','Interpreter','latex','FontSize',16);
legend(['Minimum overreach: ',num2str(Overreach.min,'%f'),'\% (step ',num2str(Overreach.minIndex),')',newline,'Maximum overreach: ',num2str(Overreach.max,'%f'),'\% (step ',num2str(Overreach.maxIndex),')'],'Interpreter','latex','location','best', 'box','off','FontSize',10)
box on
hold off

% Plot runtime for addition of new link (row 2, column 1 OR row 1, column 3 without QMinimization):
nexttile
hold on
plot(links.array,LinkRuntime.values);
if exist('QMinimization','var')
    set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(LinkRuntime.values))],'TickLabelInterpreter','latex','XAxisLocation','top');
else
    set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(LinkRuntime.values))],'TickLabelInterpreter','latex');
    xlabel('Step number','Interpreter','latex','FontSize',16);
end
ylabel('Runtime for link addition (s)','Interpreter','latex','FontSize',16);
legend(sprintf(['Minimum link runtime: ',num2str(LinkRuntime.min,'%f'),' seconds (step ',num2str(LinkRuntime.minIndex),')',newline,'Maximum link runtime: ',num2str(LinkRuntime.max,'%f'),' seconds (step ',num2str(LinkRuntime.maxIndex),')']),'Interpreter','latex','location','best', 'box','off','FontSize',10)
box on
hold off

% Plot runtime for Q_minimization (row 2, column 2 OR N/A without QMinimization):
if exist('QMinimization','var')
    nexttile
    hold on
    plot(links.array(1:length(QMinimization.values)),QMinimization.values);
    set(gca,'FontSize',10,'XLim',[0 links.tracker],'YLim',[0 (max(QMinimization.values))],'TickLabelInterpreter','latex','YAxisLocation','right','XAxisLocation','top');
    ylabel('Runtime for $Q_\mathrm{minimization}$ (s)','Interpreter','latex','FontSize',16);
    legend(['Minimum $Q_\mathrm{minimization}$ runtime: ',num2str(QMinimization.min,'%f'),' seconds (step ',num2str(QMinimization.minIndex),')',newline,'Maximum $Q_\mathrm{minimization}$ runtime: ',num2str(QMinimization.max,'%f'),' seconds (step ',num2str(QMinimization.maxIndex),')'],'Interpreter','latex','location','best', 'box','off','FontSize',10)
    box on
    hold off
end

% Format title:
sgtitle(['(',sims.objectName,'): ',sims.objectType,' discharge runtime results over ',num2str(links.tracker),' steps'],'Interpreter','latex','FontSize',24);
if exist('QMinimization','var')
    set(gcf,'Position',[0,0,800,600]);
else
    set(gcf,'Position',[0,0,1200,400]);
end
set(gcf,'Resize','off');

% Export figure:
exportgraphics(gcf,[sims.pathPNGs,'/RuntimeResults_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);