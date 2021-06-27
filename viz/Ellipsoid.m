function []=Ellipsoid(xq,yq,zq,a,b,c,color,NbRad)
% close all
% clear all
% clc
% xq = 1;
% yq = 1;
% zq = 1;
% a = 3;
% b = 2;
% c =1;
% n = 400;
[X,Y,Z] = ellipsoid(xq,yq,zq, a,b,c, NbRad);
c = ones(NbRad+1,NbRad+1);
surf(X,Y,Z,c,'FaceColor','none','EdgeColor',color);
% axis([-5 5 -5 5 -5 5])
end