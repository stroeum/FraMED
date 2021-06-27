function []=Ellipse(xq,yq,zq,ax,ay,az,color,NbVer,NbRad,NbOrthoRad)

span = floor(360/NbRad);

t    = -pi:2*pi/359:pi;

xE   = ax*cos(t);
yE   = ay*sin(t);
zE   = (zq-az/2): az/(NbVer-1) :(zq+az/2);
r    = 1/NbOrthoRad:1/NbOrthoRad:1;

hold on
for n = 1:NbOrthoRad
    plot3(xq+r(n).*xE, yq+r(n).*yE, zE(1)*ones(size(t)),'Color',color);
    plot3(xq+r(n).*xE, yq+r(n).*yE, zE(end)*ones(size(t)),'Color',color);
end

for n = 1:NbVer
    plot3(xq+xE,yq+yE,zE(n)*ones(size(t)),'Color',color);
end

for n = 1:NbRad
    plot3([xq xq+xE(n*span)],[yq yq+yE(n*span)],[zE(1) zE(1)],'Color',color);
    plot3([xq xq+xE(n*span)],[yq yq+yE(n*span)],[zE(end) zE(end)],'Color',color);
    plot3([xq+xE(n*span) xq+xE(n*span)],[yq+yE(n*span) yq+yE(n*span)],[zE(1) zE(end)],'Color',color);
end
hold off

end