% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File Name: resizeDomain.m                                               %
% Purpose: Adds capabilities to interpolate results output by TREBEC and  %
%          the respective NASA GRAM into different FraMED domain sizes    %
%          without the need to rerun any simulations. Current version     %
%          only works when the initial altitude is zero. A future version %
%          may be introduced to accommodate for alternative cases.        %
% Author: Annelisa Esparza                                                %
% Contact: annelisa.esparza@my.erau.edu                                   %
% Date Added: October 14, 2024                                              %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Select object of interest:
% Check whether a resize is being requested again or if variables have not been cleared:
if exist('grams','var')
    prompt.continue = strcat("\nDo you still wish to focus on: ",grams.objectName,' (',grams.objectType,'s)? (Y/N)\n-->');
    answer.continue = input(prompt.continue,'s');
    while ~strcmp(answer.continue,'Y') && ~strcmp(answer.continue,'N')
        fprintf('\n\tNot an acceptable input. Please enter Y for yes or N for no.\n');
        answer.continue = input(prompt.continue,'s');
    end
    if strcmp(answer.continue,'N')
        clearvars grams answer
    end
end

% Select the body of interest (and subsequently the pathway to the existing datafiles):
if ~exist('grams','var') 
    prompt1 = "\nWhat is the planetary body that the simulation is focused on?\n(No quotation marks needed for string input)\n-->";
    grams.objectName = input(prompt1,'s');
    while ~exist(grams.objectName,'dir')
        fprintf('\n\tNot an acceptable input.\n\tDefault options: Earth, Mars, Titan, Venus.\n');
        grams.objectName = input(prompt1,'s');
    end
    prompt2 = "\nWhat type of discharge should it be? (leader/streamer)\n-->";
    grams.objectType = input(prompt2,'s');
    while ~strcmp(grams.objectType,'streamer') && ~strcmp(grams.objectType,'leader')
        fprintf('\n\tNot an acceptable input. Please enter streamer or leader.\n');
        grams.objectType = input(prompt2,'s');
    end
    grams.ogObjectType = "leader";
end 

% Loads current altitude values for the respective object:
grams.z = load(strcat(grams.objectName,'/',grams.objectName,'_z.dat'));
grams.Nz = size(grams.z,1);
grams.Lmin = grams.z(1);
if grams.Lmin ~= 0
    fprintf('The resizeDomain.m script is currently not valid for nonzero starting altitudes.\n');
    return;
end
grams.Lmax = grams.z(grams.Nz);

% Outputs relevant information about the existing simulation domain's parameters:
spacingCheck = norm((grams.Lmin:(grams.z(2)-grams.z(1)):grams.Lmax)'-grams.z);
if spacingCheck==0
    grams.dz = grams.z(2)-grams.z(1);
    fprintf(['\n*** Domain for ',grams.objectName,' has an altitude domain between ',num2str(grams.Lmin), 'm --> ',num2str(grams.Lmax),'m with ',num2str(grams.dz),'m spacings in between (total of ',num2str(grams.Nz),' grid points). ***\n']);
else
    fprintf(['\n*** Domain for ',grams.objectName,' has an altitude domain between ',num2str(grams.Lmin), 'm --> ',num2str(grams.Lmax),'m with uneven spacings in between (total of ',num2str(grams.Nz),' grid points). ***\n']);
end

%% Determine whether to change spacings, domain size, or both:
prompt.changes = "\nWould you like to change the spacings (D), maximum domain altitude (L), or both (B)?\n-->";
answer.changes = input(prompt.changes,'s');
while ~strcmp(answer.changes,'D') && ~strcmp(answer.changes,'L') && ~strcmp(answer.changes,'B') && (strcmp(answer.changes,'L') && spacingCheck~=0)
    if strcmp(answer.changes,'L') && spacingCheck~=0
        fprintf('\n\tCannot define new domain without even spacings. Enter B to change the domain size and define the spacings.\n');
    else
        fprintf('\n\tNot an acceptable input. Please enter D for spacings, L for maximum altitude, or B to change both.\n');
    end
    answer.changes = input(prompt.changes,'s');
end

%% Clarify specificities of request:
switch answer.changes
    case 'D'
        fprintf('\n*** (D) CHANGING SPACING BETWEEN ALTITUDES ***\n');
        % Requests user-input for new spacings distance:
        if spacingCheck==0
            prompt.spacings = strcat("\nInput the new spacing value in meters (currently ",num2str(grams.dz)," meters).\n-->");
        else
            prompt.spacings = strcat("\nInput the new spacing value in meters (currently non-uniformly spaced).\n-->");
        end
        answer.spacings = input(prompt.spacings);
        answer.domain = grams.Lmax;

        % Ensures spacing and altitude values align:
        while mod(answer.domain,answer.spacings)~=0
            fprintf(strcat("\n\tMaximum of ",num2str(answer.domain),"m cannot be achieved with requested spacings.\n"));
            if answer.domain >= grams.Lmax
                suggestedmax = floor(grams.Lmax/answer.spacings)*answer.spacings;
            else
                suggestedmax = round(answer.domain/answer.spacings)*answer.spacings;
            end
            prompt.confirm = strcat("\nIs a maximum altitude of ",num2str(suggestedmax)," meters (multiple of ",num2str(answer.spacings),"m spacings) acceptable? (Y/N)\n-->");
            answer.confirm = input(prompt.confirm,'s');
            while ~strcmp(answer.confirm,'N') && ~strcmp(answer.confirm,'Y')
                fprintf('\n\tNot an acceptable input. Please enter Y for yes or N for no.\n');
                answer.confirm = input(prompt.confirm,'s');
            end
            if strcmp(answer.confirm,'Y')
                answer.domain = suggestedmax;
            elseif strcmp(answer.confirm,'N')
                prompt.domain = strcat("\nInput the new maximum altitude in meters (currently incompatible with ",num2str(answer.domain)," meters).\n-->");
                answer.domain = input(prompt.domain);
            end
        end
    case 'L'
        fprintf('\n*** (L) CHANGING MAXIMUM ALTITUDE ***\n');
        % Requests user-input for new domain maximum:
        prompt.domain = strcat("\nInput the new maximum altitude in meters (currently ",num2str(grams.Lmax)," meters).\n-->");
        answer.domain = input(prompt.domain);
        answer.spacings = grams.dz;

        % Ensures spacing and altitude values align:
        while mod(answer.domain,answer.spacings)~=0 || answer.domain>grams.Lmax
            if answer.domain>grams.Lmax 
                fprintf('\n\tCannot extrapolate values for a larger domain size.\n');
            end
            if mod(answer.domain,answer.spacings)~=0
                fprintf('\n\tRequested maximum cannot be achieved with current spacings.\n');
                if (round(answer.domain/answer.spacings)*answer.spacings)>=grams.Lmax
                    suggestedmax = grams.Lmax;
                else
                    suggestedmax = round(answer.domain/answer.spacings)*answer.spacings;
                end
            end
            prompt.confirm = strcat("\nIs a maximum altitude of ",num2str(suggestedmax)," meters (multiple of ",num2str(answer.spacings),"m spacings) acceptable? (Y/N)\n-->");
            answer.confirm = input(prompt.confirm,'s');
            while ~strcmp(answer.confirm,'N') && ~strcmp(answer.confirm,'Y')
                fprintf('\n\tNot an acceptable input. Please enter Y for yes or N for no.\n');
                answer.confirm = input(prompt.confirm,'s');
            end
            if strcmp(answer.confirm,'Y')
                answer.domain = suggestedmax;
            elseif strcmp(answer.confirm,'N')
                prompt.domain = strcat("\nInput the new maximum altitude in meters (currently incompatible with ",num2str(answer.domain)," meters).\n-->");
                answer.domain = input(prompt.domain);
            end
        end
    case 'B'
        fprintf('\n*** (B) CHANGING MAXIMUM ALTITUDE AND SPACINGS IN-BETWEEN ***\n');
        % Requests user-input for new spacings distance:
        if spacingCheck==0
            prompt.spacings = strcat("\nInput the new spacing value in meters (currently ",num2str(grams.dz)," meters).\n-->");
        else
            prompt.spacings = strcat("\nInput the new spacing value in meters (currently non-uniformly spaced).\n-->");
        end
        answer.spacings = input(prompt.spacings);

        % Requests user-input for new domain maximum:
        prompt.domain = strcat("\nInput the new maximum altitude in meters (currently ",num2str(grams.Lmax)," meters).\n-->");
        answer.domain = input(prompt.domain);

        % Ensures spacing and altitude values align:
        while mod(answer.domain,answer.spacings)~=0
            fprintf(strcat("\n\tMaximum of ",num2str(answer.domain),"m cannot be achieved with requested spacings.\n"));             
            if (round(answer.domain/answer.spacings)*answer.spacings)>=(grams.Lmax)
                suggestedmax = floor(grams.Lmax/answer.spacings)*answer.spacings;
            else
                suggestedmax = round(answer.domain/answer.spacings)*answer.spacings;
            end
            prompt.confirm = strcat("\nIs a maximum altitude of ",num2str(suggestedmax)," meters (multiple of ",num2str(answer.spacings),"m spacings) acceptable? (Y/N)\n-->");
            answer.confirm = input(prompt.confirm,'s');
            while ~strcmp(answer.confirm,'N') && ~strcmp(answer.confirm,'Y')
                fprintf('\n\tNot an acceptable input. Please enter Y for yes or N for no.\n');
                answer.confirm = input(prompt.confirm,'s');
            end
            if strcmp(answer.confirm,'Y')
                answer.domain = suggestedmax;
            elseif strcmp(answer.confirm,'N')
                prompt.domain = strcat("\nInput the new maximum altitude in meters (currently incompatible with ",num2str(answer.domain)," meters).\n-->");
                answer.domain = input(prompt.domain);
            end  
        end    
end
%% Adapt values per request:
if strcmp(grams.objectType,'streamer') && strcmp(grams.ogObjectType,'leader')
    negpropfactor = (1/3)/(1/20);
    pospropfactor = (1/6)/(1/20);
    initfactor = 10;
elseif strcmp(grams.objectType,'leader') && strcmp(grams.ogObjectType,'streamer')
    negpropfactor = (1/20)/(1/3);
    pospropfactor = (1/20)/(1/6);
    initfactor = (1/10);
elseif strcmp(grams.objectType,grams.ogObjectType)
    negpropfactor = 1;
    pospropfactor = 1;
    initfactor = 1;
end

%% Interpolate and export:
% Summarize new altitude domain & export results:
newvals_z = (grams.Lmin:answer.spacings:answer.domain)';
newvals_N = [0; 0; size(newvals_z,1)];
fprintf(['\n*** New domain for ',grams.objectName,' has an altitude domain between ',num2str(newvals_z(1)), 'm --> ',num2str(answer.domain),'m with ',num2str(answer.spacings),'m spacings in between (total of ',num2str(newvals_N(3)),' grid points). ***\n']);

% Creates descriptive filetag prefix:
grams.pathOutput = [grams.objectName,'/',grams.objectType,'-N',num2str(size(newvals_z,1)),'-D',num2str(answer.spacings),'m'];
%grams.pathOutput = strcat(grams.pathOutput,'_half');
save([grams.pathOutput,'_z.dat'],'newvals_z','-ascii');

% Converts & saves particle density (ng), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_ng.dat'),"file") 
    grams.ng = load(strcat(grams.objectName,'/',grams.objectName,'_ng.dat'));
    newvals_ng = interp1(grams.z,grams.ng,newvals_z,'spline');
    save([grams.pathOutput,'_ng.dat'],'newvals_ng','-ascii');
else
    fprintf(strcat("\tParticle density file not found for ",grams.objectName,".\n"));
end

% Converts & saves positive propagation threshold (Eth+), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_Eth_positive_Vm.dat'),"file") 
    grams.Ethpos = pospropfactor.*load(strcat(grams.objectName,'/',grams.objectName,'_Eth_positive_Vm.dat'));
    newvals_Ethpos = interp1(grams.z,grams.Ethpos,newvals_z,'spline');
    save([grams.pathOutput,'_Eth_positive_Vm.dat'],'newvals_Ethpos','-ascii');
else
    fprintf(strcat("\tPositive propagation threshold file not found for ",grams.objectName,".\n"));
end

% Converts & saves negative propagation threshold (Eth-), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_Eth_negative_Vm.dat'),"file") 
    grams.Ethneg = negpropfactor.*load(strcat(grams.objectName,'/',grams.objectName,'_Eth_negative_Vm.dat'));
    newvals_Ethneg = interp1(grams.z,grams.Ethneg,newvals_z,'spline');
    save([grams.pathOutput,'_Eth_negative_Vm.dat'],'newvals_Ethneg','-ascii');
else
    fprintf(strcat("\tNegative propagation threshold file not found for ",grams.objectName,".\n"));
end

% Converts & saves initiation threshold (E_initiation), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_E_initiation_Vm.dat'),"file") 
    grams.Einit = initfactor.*load(strcat(grams.objectName,'/',grams.objectName,'_E_initiation_Vm.dat'));
    newvals_Einit = interp1(grams.z,grams.Einit,newvals_z,'spline');
    save([grams.pathOutput,'_E_initiation_Vm.dat'],'newvals_Einit','-ascii');
else
    fprintf(strcat("\tInitiation threshold file not found for ",grams.objectName,".\n"));
end

% Converts & saves temperature (Tg), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_Tg.dat'),"file") 
    grams.Tg = load(strcat(grams.objectName,'/',grams.objectName,'_Tg.dat'));
    newvals_Tg = interp1(grams.z,grams.Tg,newvals_z,'spline');
    save([grams.pathOutput,'_Tg.dat'],'newvals_Tg','-ascii');
else
    fprintf(strcat("\tTemperature file not found for ",grams.objectName,".\n"));
end

% Converts & saves reduced electric field threshold (Ek), if the file exists:
if exist(strcat(grams.objectName,'/',grams.objectName,'_Ek.dat'),"file") 
    grams.Ek = load(strcat(grams.objectName,'/',grams.objectName,'_Ek.dat'));
    newvals_Ek = interp1(grams.z,grams.Ek,newvals_z,'spline');
    save([grams.pathOutput,'_Ek.dat'],'newvals_Ek','-ascii');
else
    fprintf(strcat("\tReduced electric field file not found for ",grams.objectName,".\n"));
end

%% (Optional) Define total (xyz) domain size for read-in:
prompt.domdef = strcat("\nWould you like to define the number of grid points for the x and y dimensions too?\nThis is not required but simplifies the use of main.cpp to initialize the simulation. (Y/N)\n-->");
answer.domdef = input(prompt.domdef,'s');
while ~strcmp(answer.domdef,'Y') && ~strcmp(answer.domdef,'N')
    fprintf('\n\tNot an acceptable input. Please enter Y for yes or N for no.\n');
    answer.domdef = input(prompt.domdef,'s');
end

% Inquires and saves spacing and grid values:
if strcmp(answer.domdef,'Y')
    % Determine x-dimension sizing:
    prompt.xnodes = strcat("\nThere are currently ",num2str(newvals_N(3))," nodes for the z-dimension.\nHow many grid points shall the x-dimension span?\n-->");
    newvals_N(1) = input(prompt.xnodes);
    while newvals_N(1)~=round(newvals_N(1))
        fprintf('\n\tValue must be an integer.\n');
        newvals_N(1) = input(prompt.xnodes);
    end

    % Determine y-dimension sizing:
    prompt.ynodes = strcat("\nHow many grid points shall the y-dimension span?\n-->");
    newvals_N(2) = input(prompt.ynodes);
    while newvals_N(2)~=round(newvals_N(2))
        fprintf('\n\tValue must be an integer.\n');
        newvals_N(2) = input(prompt.ynodes);
    end

    % Defines total domain size with assumption of equal spacings:
    newvals_D = [answer.spacings; answer.spacings; answer.spacings];
    newvals_L = newvals_D.*(newvals_N-1);
   
    % Outputs summary to the screen:
    fprintf(strcat("\n\t******* SUMMARY OF INITIALIZED SIMULATION DOMAIN *******" + ...
        "\n\t\t\t\t(x)\t(y)\t(z)" + ...
        "\n\tDomain Size \t(L):\t",num2str(newvals_L(1)),"\t",num2str(newvals_L(2)),"\t",num2str(newvals_L(3)),"\t(meters)" + ...
        "\n\tSpacings \t(D):\t",num2str(newvals_D(1)),"\t",num2str(newvals_D(2)),"\t",num2str(newvals_D(3)),"\t(meters)" + ...
        "\n\tGrid Points \t(N):\t",num2str(newvals_N(1)),"\t",num2str(newvals_N(2)),"\t",num2str(newvals_N(3)),"\t(nodes)\n"));
    if exist(strcat(grams.objectName,'/',grams.objectName,'_Nxyz.dat'),"file") 
        old_N = load(strcat(grams.objectName,'/',grams.objectName,'_Nxyz.dat'));
        reduction = (newvals_N(1)*newvals_N(2)*newvals_N(3))/(old_N(1)*old_N(2)*old_N(3));
        fprintf(strcat("\n\tTotal nodes = ",num2str(newvals_N(1)*newvals_N(2)*newvals_N(3),'%d'),"\t(",num2str(100*reduction,'%.1f'),"%% of previous ",num2str(old_N(1)*old_N(2)*old_N(3),'%d'),", should run ",num2str(round(10/reduction)/10,'%.1f'),"x faster)\n"));
    else
        fprintf(strcat("\n\tTotal nodes = ",num2str(newvals_N(1)*newvals_N(2)*newvals_N(3),'%1.2e'),"\n"));
    end

    % Outputs the path and prefix for the output files:
    fprintf(strcat("\n\tPath and prefix of associated files: ",grams.pathOutput,"_\n"));

    % Exports grid size and spacings:
    save([grams.pathOutput,'_Nxyz.dat'],'newvals_N','-ascii');
    save([grams.pathOutput,'_Dxyz.dat'],'newvals_D','-ascii');
else
    % Outputs the path and prefix for the output files:
    fprintf(strcat("\n\tPath and prefix of associated files: ",grams.pathOutput,"_\n"));
end