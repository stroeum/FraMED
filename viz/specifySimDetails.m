function sims = specifySimDetails()
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
    sims.pathEPSs = ['../Figures/',sims.objectName,'/',sims.objectType,'/EPSs'];
    if ~exist(sims.pathEPSs,'dir')
        mkdir(sims.pathEPSs);
    end

    % Specifies the boundary conditions for the simulation:
    prompt_BCtype = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
    sims.BCtype = input(prompt_BCtype,'s');                    
    while ~strcmp(sims.BCtype,'FS') && ~strcmp(sims.BCtype,'G')
        fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
        sims.BCtype = input(prompt_BCtype,'s');
    end
end