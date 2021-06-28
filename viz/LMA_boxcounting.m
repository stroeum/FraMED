%% Init
clear all
close all
clc

%% xy
c = imread('xy.png');

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

%% xy
c = imread('xz.png');

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
c = imread('yz.png');

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

print(gcf,'-depsc','fractal_box')