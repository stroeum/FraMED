%% Init
clearvars -except sims
close all


if ~exist('sims','var') || ~isfield(sims,'pathPNGs') || ~isfield(sims,'pathVideos')
    specifySimDetails;
end

if ~exist(sim.pathPNGs, 'dir')
    fprintf('Files required do not exist.');
    return;
end

%% xy
c = imread(strcat(sims.pathPNGs,'/xy_',sims.objectName,'_',sims.objectType,'.png'));

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
c = imread(strcat(sims.pathPNGs,'/xz_',sims.objectName,'_',sims.objectType,'.png'));

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
c = imread(strcat(sims.pathPNGs,'/yz_',sims.objectName,'_',sims.objectType,'.png'));

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