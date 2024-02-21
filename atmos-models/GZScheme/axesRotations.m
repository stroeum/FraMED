function [xRotation, yRotation] = axesRotations(azDeg,elDeg)
    projectedunitx = rotx(elDeg) * rotz(-azDeg) * [1;0;0];
    projectedunity = rotx(elDeg) * rotz(-azDeg) * [0;1;0];
    xRotation = atan2d(projectedunitx(3),projectedunitx(1));
    yRotation = atan2d(projectedunity(3),projectedunity(1));
end

function Rx = rotx(angle)
    Rx = [1 0 0; 0 cosd(angle) -sind(angle); 0 sind(angle) cosd(angle)];
end

function Rz = rotz(angle)
    Rz = [cosd(angle) -sind(angle) 0; sind(angle) cosd(angle) 0; 0 0 1];
end