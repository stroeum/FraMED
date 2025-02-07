function sims = specifySimDetails()
    prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
    sims.objectName = input(prompt1,'s');

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

    % Specifies the boundary conditions for the simulation:
    if exist('../results/Type_Discharge.txt','file')
        sims.BCtype = string(textscan(fopen('../results/Type_BC.txt'),'%s'));
    else        
        prompt_BCtype = '\nIs the domain in free space (FS) or is z = 0 grounded (G)?\n-->';
        sims.BCtype = input(prompt_BCtype,'s');                    
        while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
            fprintf('\n\tNot an acceptable input. Please enter FS (for free space) or G (for grounded).\n');
            sims.BCtype = input(prompt_BCtype,'s');
        end
    end

    if exist('../results/Type_Result.txt','file')
        sims.disType = string(textscan(fopen('../results/Type_Result.txt'),'%s'));
    end

    % Settings to ensure proper directory referencing:
    sims.pathPNGs = strcat('../Figures/',sims.objectName,'/',sims.objectType,'/PNGs')
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
end