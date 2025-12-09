% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  File Name: Plot3D_CloudDistribution.m                                  %
%    Purpose: Visualizes the charged cloud structure and ground with      %
%             custom colormap. Outputs a figure to the screen but does    %
%             not save it to a file automatically.                        %
%     Author: Annelisa Esparza                                            %
%    Contact: annelisa.esparza@my.erau.edu                                %
% Added Date: February 22, 2022                                           %
%    Updates: February 2024 - Integrated free-space BCs.                  %
%              October 2025 - Integrated checkMagnitude.m function and    %
%                             tin-can BCs.                                %
%             December 2025 - Integrated the changes made to the          %
%                             specifySimDetails.m function and the newly  %
%                             introduced setUpAxes.m function.            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

close all
clearvars -except sims

%% User-defined variables that influence the plot:
saveAfter = 'N';

%% Fully-automated onwards:
clf
fprintf('\n*** Executing Plot3D_CloudDistribution.m script. ***\n');

if ~exist('sims','var')
    specifySimDetails;
end

% Load data files
cd ../results
if strcmp(saveAfter,'Y')
    fprintf('\tUsing charge layer data from after the discharge occurs.\n');
    rho.data = load('rho3d.dat,'-ascii'); % After discharge
else
    fprintf('\tUsing charge layer data from before the discharge occurs.\n');
    rho.data = load('rhoAmb.dat','-ascii'); % Before discharge
end
cd ../viz
rho.data = convertTo3d(rho.data,sims); % _nC/_m^3

% Linear spaces for the three position dimensions:
x = (0:sims.domain.dx:sims.domain.maxx)'*sims.spatialFactor.Number;
y = (0:sims.domain.dy:sims.domain.maxy)'*sims.spatialFactor.Number;
z = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)'*sims.spatialFactor.Number;

[X,Y,Z] = meshgrid(x,y,z);
rho.max = max(max(max(rho.data)));
rho.min = min(min(min(rho.data)));

% Plotting figure:
figure(1);
set(gcf,'Position',[0,0,sims.plotWidth,sims.plotHeight]);
hold on;

% Calls the function that automatically recognizes charge regions:
plottingChargeRegions('scalar',0.4,rho,X,Y,Z,sims);

% Save figure as file:
if strcmp(saveAfter,'Y')
    exportgraphics(gcf,strcat(sims.pathPNGs,'/CloudDistribution_After_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);  % Before discharge
else
    exportgraphics(gcf,strcat(sims.pathPNGs,'/CloudDistribution_Before_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);   % After discharge
end