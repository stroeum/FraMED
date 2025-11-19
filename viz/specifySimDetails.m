% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: specifySimDetails.m                                         %
%    Purpose: Prompts user to define the object simulated, but otherwise  %
%             automatically determines the type of discharge simulated    %
%             based on files in the 'results' directory. Resulting 'sims' %
%             structure then informs the output directory for any figures %
%             as well as certain aspects of the figures (e.g. any         %
%             grounding regions based on boundary conditions, figure      %
%             titles, propagation speeds, etc).                           %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: December 1, 2024                                            %
%    Updates: February 2025 - Automated the definition of leader/streamer %
%                             discharge types.                            %
%                 June 2025 - Integrated the automated read-in of         %
%                             simulation conditions following the         % 
%                             introduction of output txt files from the   %
%                             C++ repository.                             %
%              October 2025 - Expanded the boundary condition options.    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if ~exist('sims','var') || (exist('sims','var') && ~isfield(sims,'objectName'))
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
end

% Determines whether the simulation is for a leader or streamer:
if exist('../results/Type_Discharge.txt','file')
    sims.objectType = string(textscan(fopen('../results/Type_Discharge.txt'),'%s'));
else
    prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
    sims.objectType = input(prompt2,'s');
    while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
        fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
        sims.objectType = input(prompt2,'s');
    end
end

% Consolidate the definitions for the associated propagation speeds:
if strcmp(sims.objectType, 'Leader')
    sims.vprop.pos = 4.4*(10^5); % propagation speed for positive leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
    sims.vprop.neg = 4.4*(10^5); % propagation speed for negative leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
elseif strcmp(sims.objectType, 'Streamer')
    sims.vprop.pos = 3.0*(10^7); % propagation speed for positive streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9
    sims.vprop.neg = 3.0*(10^7); % propagation speed for negative streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9
end

% Specifies the boundary conditions for the simulation:
if exist('../results/Type_Discharge.txt','file')
    sims.BCtype = string(textscan(fopen('../results/Type_BC.txt'),'%s'));
else        
    prompt_BCtype = '\nIs the domain in free space (FS), is z = 0 grounded (G), are the top and bottom of the domain grounded (G_G), or are all sides of the domain grounded (TIN_CAN)?\n-->';
    sims.BCtype = input(prompt_BCtype,'s');                    
    while ~strcmp(sims.BCtype,'FS') && ~strcmp(sims.BCtype,'G') && ~strcmp(sims.BCtype,'TIN_CAN') && ~strcmp(sims.BCtype,'G_G')
        fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
        sims.BCtype = input(prompt_BCtype,'s');
    end
end

if exist('../results/Type_Result.txt','file')
    sims.disType = string(textscan(fopen('../results/Type_Result.txt'),'%s'));
end

% Settings to ensure proper directory referencing:
sims.pathPNGs = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/PNGs');
if ~exist(sims.pathPNGs,'dir')
    mkdir(sims.pathPNGs);
end
sims.pathVideos = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/Videos');
if ~exist(sims.pathVideos,'dir')
    mkdir(sims.pathVideos);
end
sims.pathEPSs = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/EPSs');
if ~exist(sims.pathEPSs,'dir')
    mkdir(sims.pathEPSs);
end