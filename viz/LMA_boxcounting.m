%% Init
clearvars -except sims
close all


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
end 
exportgraphics(gcf,[sims.pathPNGs,'/FieldLines_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);

if ~exist(sim.pathPNGs, 'dir')
    fprintf('Files required do not exist.');
    return;
end

%% xy
c = imread([sims.pathPNGs,'/xy_',sims.objectName,'_',sims.objectType,'.png']);

figure(1)
subplot(331)
imagesc(~c);
colormap gray
axis image

subplot(332)
boxcount(c,'slope')

[n, r] = boxcount(c);
subplot(333)
loglog(r, n,'bo-', r, (r/r(end)).^(-2), 'r--')
xlabel('r')
ylabel('n(r)')
legend('actual box-count','space-filling box-count');

%% xz
c = imread([sims.pathPNGs,'/xz_',sims.objectName,'_',sims.objectType,'.png']);

subplot(334)
imagesc(~c);
colormap gray
axis image

subplot(335)
boxcount(c,'slope')

[n, r] = boxcount(c);
subplot(336)
loglog(r, n,'bo-', r, (r/r(end)).^(-2), 'r--')
xlabel('r')
ylabel('n(r)')
legend('actual box-count','space-filling box-count');

%% yz
c = imread([sims.pathPNGs,'/yz_',sims.objectName,'_',sims.objectType,'.png']);

subplot(337)
imagesc(~c);
colormap gray
axis image

subplot(338)
boxcount(c,'slope')

[n, r] = boxcount(c);
subplot(339)
loglog(r, n,'bo-', r, (r/r(end)).^(-2), 'r--')
xlabel('r')
ylabel('n(r)')
legend('actual box-count','space-filling box-count');

print(gcf,'-depsc','../Figures/fractal_box')