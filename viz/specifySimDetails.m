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
%             November 2025 - Included propagation speed definitions.     %
%             December 2025 - Integrated recurrent domain information.    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Determines whether the simulation is for a leader or streamer:
if exist('../results/Type_Discharge.txt','file') && exist('../results/Type_BC.txt','file')
    sims.objectType = string(textscan(fopen('../results/Type_Discharge.txt'),'%s'));
    sims.BCtype     = string(textscan(fopen('../results/Type_BC.txt'),'%s'));
else
    error('Simulation is undefined. Ensure there are results to pull from.');
end

% Consolidates the definitions for the associated propagation speeds:
if strcmp(sims.objectType, 'Leader')
    sims.vprop.pos = 4.4*(10^5); % propagation speed for positive leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
    sims.vprop.neg = 4.4*(10^5); % propagation speed for negative leaders,   ~440 km/s,    Thomson1985, doi:10.1029/JD090iD05p08136
elseif strcmp(sims.objectType, 'Streamer')
    sims.vprop.pos = 3.0*(10^7); % propagation speed for positive streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9
    sims.vprop.neg = 3.0*(10^7); % propagation speed for negative streamers, ~30,000 km/s, Pasko2012,   doi:10.1007/s11214-011-9813-9
end

% Preloads and stores attributes regarding the domain:
if exist('../results/dxyz.dat','file') && exist('../results/Nxyz.dat','file') 
    load('../results/Nxyz.dat');
    load('../results/dxyz.dat');
    if exist('../results/z_gnd.dat','file')
        sims.domain.gnd = load('../results/z_gnd.dat');
    else
        sims.domain.gnd = 0;
    end
    sims.domain.Nx = Nxyz(1);
    sims.domain.Ny = Nxyz(2);
    sims.domain.Nz = Nxyz(3);
    sims.domain.dx = dxyz(1);
    sims.domain.dy = dxyz(2);
    sims.domain.dz = dxyz(3);
    clear Nxyz dxyz
    sims.domain.minx = 0;
    sims.domain.maxx = (sims.domain.Nx-1)*sims.domain.dx;
    sims.domain.miny = 0;
    sims.domain.maxy = (sims.domain.Ny-1)*sims.domain.dy;
    sims.domain.minz = 0;
    sims.domain.maxz = ((sims.domain.Nz-1)*sims.domain.dz)+sims.domain.gnd;

    % Determines the spatial factor for the domain:
    sims.spatialFactor = checkMagnitude([(0:sims.domain.dy:sims.domain.maxx)'; (0:sims.domain.dy:sims.domain.maxy)'; (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)']);

    % Assigns plot height and width based on domain:
    sims.plotWidth = 600;
    sims.plotHeight = round((sims.domain.Nz/max([sims.domain.Nx sims.domain.Ny]))*5)*80;
else
    error('Simulation domain is undefined. Ensure there are results to pull from.');
end

% Specifies the status of the discharge:
if exist('../results/Type_Result.txt','file')
    sims.disType = string(textscan(fopen('../results/Type_Result.txt'),'%s'));
else
    sims.disType = "Ongoing";
end

% Inquires the simulation object (determines naming convention):
if ~exist('sims','var') || (exist('sims','var') && ~isfield(sims,'objectName'))
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');
end

% Settings to ensure proper directory referencing:
baseDirectory = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/',sims.disType);
copyfile('../src/main.cpp',baseDirectory);
fprintf(strcat("\nVisualizations for this simulation will be stored within the following directory --> ",baseDirectory,"\n\n"));
sims.pathPNGs = strcat(baseDirectory,'/PNGs');
if ~exist(sims.pathPNGs,'dir')
    mkdir(sims.pathPNGs);
end
sims.pathVideos = strcat(baseDirectory,'/Videos');
if ~exist(sims.pathVideos,'dir')
    mkdir(sims.pathVideos);
end
sims.pathEPSs = strcat(baseDirectory,'/EPSs');
if ~exist(sims.pathEPSs,'dir')
    mkdir(sims.pathEPSs);
end