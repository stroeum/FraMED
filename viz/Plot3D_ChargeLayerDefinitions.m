% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   File Name: Plot3D_ChargeLayerDefinitions.m                            %
%     Purpose: Visualizes the charged cloud structure and ground with     %
%              custom colormap. Labels the charge layers based on the     %
%              inputs for FraMED.                                         %
%      Author: Annelisa Esparza                                           %
%     Contact: annelisa.esparza@my.erau.edu                               %
%  Added Date: April 29, 2022                                             %
%     Updates: February 2025 - Updated to match the options available for %
%                              the createCustomColorMap function.         %  
%               October 2025 - Integrates the checkMagnitude.m function.  %
%              December 2025 - Integrated the changes made to the         %
%                              specifySimDetails.m function and the newly %
%                              introduced setUpAxes.m function.           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Initiate
%close all
clearvars -except sims
figure(1)
clf
fprintf('\n*** Executing Plot3D_ChargeLayerDefinitions.m script. ***\n');
if ~exist('sims','var')
    specifySimDetails;
end

%% Load data files
cd ../results
rho.data = load('rhoAmb.dat','-ascii');
load('ChargeLayers.dat',     '-ascii');
cd ../viz

rho.data = convertTo3d(rho.data,sims); % _nC/_m^3

clear InitPoint

Q = ChargeLayers(:,1);
R = ChargeLayers(:,5)*sims.spatialFactor.Number;
h = ChargeLayers(:,7)*sims.spatialFactor.Number;
center = ChargeLayers(:,2:4)*sims.spatialFactor.Number;

% Linear spaces for the three position dimensions:
x = (0:sims.domain.dx:sims.domain.maxx)'*sims.spatialFactor.Number;
y = (0:sims.domain.dy:sims.domain.maxy)'*sims.spatialFactor.Number;
z = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)'*sims.spatialFactor.Number;

[X,Y,Z] = meshgrid(x,y,z);
rho.max = max(max(max(rho.data)));
rho.min = min(min(min(rho.data)));

%% Plotting figure:
% Draw the tree
figure(1);
set(gcf,'Position',[0,0,1.5*sims.plotWidth,1.5*sims.plotHeight]);
hold on;

% Calls the function that automatically recognizes charge regions:
plottingLayerDefs('scalar',0.4,rho,X,Y,Z,R,h,Q,center,sims);
setUpAxes(sims,'xyz');
legend
title(strcat("Charge Layer Distribution (",sims.disType," ",sims.objectType," on ",sims.objectName,")"),'Interpreter','latex','FontSize',28);
exportgraphics(gcf,strcat(sims.pathPNGs,'/ChargeLayerDefs_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);

