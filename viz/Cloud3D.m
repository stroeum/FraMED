function []=Cloud3D(zB,zT, Lx, Ly, Lz, z_gnd)
% z.T = 15;
% z.B = 3;
% N.z = 23;
z.B = zB+z_gnd;
z.T = zT+z_gnd;
N.z = 23;
d.z = (z.T-z.B)/(N.z-1);

z.o = z.B:d.z:z.T;
R.o = [0 5 3 4.5 3.5 3.25 6 5.5 5.75 4.5 5 4.5 5.5 4.5 5.75 6.25 5.5 4 7 9 7 3.5 0]*7.5/10*1e3;

z.i = z.B:d.z/100:z.T;
R.i = interp1(z.o,R.o,z.i,'spline');
x             = R.i;
y             = z.i;
[theta,rho,z] = cart2pol(x, zeros(size(x,1), size(x,2)), y);
r             = 25;
theta         = repmat( linspace(0,2*pi,r)', 1, length(rho(:)) );
rho           = repmat(rho(:).', r, 1);
z             = repmat(z(:).', r, 1);
[X,Y,Z]       = pol2cart(theta,rho,z);

C = ones(size(Z,1),size(Z,2));
surf((X+Lx/2)*1e-3,(Y+Ly/2)*1e-3,(Z)*1e-3,C,'EdgeColor','none','FaceColor',[.75 .75 .75],'FaceAlpha',.5);
axis([0 Lx 0 Ly 0 (Lz+z_gnd)]*1e-3)
axis equal
view(3)

end