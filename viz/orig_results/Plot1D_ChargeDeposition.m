% Lightning charge deposition as a function of altitude:
% Analysis of the charge neutrality shows slight differences as output of
% C++ and Matlab. In particular the positive and negative channel charges
% induced on the channel, may be slightly different. This is mostly due to
% the implicit conversion from double to float when data derived from the
% C++ simulation to text files. As in the C++ code, we only account for the
% charge at the channel location. If this is not done, noise exists in
% rho_cha, in particular due to the boundaries update which introduces
% artificial gradients of potential and hence charge density close to the
% boundaries.

close all
clear all
clc
beep  off

%% Load relevant data
dxyz              = load('dxyz.dat')*1e-3;                                 %_km
Nxyz              = load('Nxyz.dat');                                      %_dimensionless
Links.data        = load('EstablishedLinks.dat');
z_gnd             = load('z_gnd.dat')*1e-3;                                %_km
ChargeLayers.data = load('ChargeLayers.dat')*1e-3;                         %_kC, _km
rho.tot           = load('rho.dat');                                       %_nC/m^3
rho.amb           = load('rhoAmb.dat');                                    %_nC/m^3

%% Convert to 3D
rho.amb           = ConvertTo3d(rho.amb,Nxyz);
rho.tot           = ConvertTo3d(rho.tot,Nxyz);
rho.cha           = rho.tot-rho.amb;

%% Derive main parameters
Nb.Links          = size(Links.data);
Nb.Links          = Nb.Links(1);
Nb.Points         = Nb.Links+1;
Nb.ChargeLayers   = size(ChargeLayers.data);
Nb.ChargeLayers   = Nb.ChargeLayers(1);

N.x               = Nxyz(1);
N.y               = Nxyz(2);
N.z               = Nxyz(3);

d.x               = dxyz(1);                                               %_km
d.y               = dxyz(2);                                               %_km
d.z               = dxyz(3);                                               %_km

L.x               = (N.x-1)*d.x;                                           %_km
L.y               = (N.y-1)*d.y;                                           %_km
L.z               = (N.z-1)*d.z;                                           %_km
clear Nxyz dxyz InitPoint

ChargeLayers.Type       = 'disks'; % 'spheres'; %
ChargeLayers.Line.Style = '--';
ChargeLayers.Line.Width = 1;
ChargeLayers.Edge.Color = [[.75 .75 .75];[.75 .75 .75];[.75 .75 .75];[.75 .75 .75]];
% ChargeLayers.Edge.Color = ['none';'none';'none';'none'];
Ground.Line.Width       = .25;

%% Charge deposition
% Q.cha.plus(N.z)  = 0;
% Q.cha.minus(N.z) = 0;
% for kk=2:N.z-1
%     for jj=2:N.y-1
%         for ii=2:N.x-1
%             if(rho.cha(ii,jj,kk)>=0)
%                 Q.cha.plus(kk)  = Q.cha.plus(kk)  + rho.cha(ii,jj,kk)*d.x*d.y*d.z;
%             elseif(rho.cha(ii,jj,kk)<=0)
%                 Q.cha.minus(kk) = Q.cha.minus(kk) + rho.cha(ii,jj,kk)*d.x*d.y*d.z;
%             end
%         end
%     end
% end

%% Zoom Area
FocusArea.x(1)    = 0;
FocusArea.x(2)    = L.x;
FocusArea.y(1)    = 0;
FocusArea.y(2)    = L.y;
FocusArea.z(1)    = z_gnd;
FocusArea.z(2)    = L.z+z_gnd;

%% Fonts
Font.Size.Labels  = 12;
Font.Size.Axis    = 10;
Font.Name         = 'Helvetica';

%% Altitude histogram
alt_histogram.data(N.z)               = 0;
Q.cha.plus(N.z)                       = 0;
Q.cha.minus(N.z)                      = 0;

alt_histogram.data(Links.data(1,3)+1) = 1;
d.q = rho.cha(Links.data(1,1)+1,Links.data(1,2)+1,Links.data(1,3)+1) *d.x*d.y*d.z;
if (d.q>=0)
    Q.cha.plus(Links.data(1,3)+1)  = Q.cha.plus(Links.data(1,3)+1)  + d.q;
elseif (d.q<=0)
    Q.cha.minus(Links.data(1,3)+1) = Q.cha.minus(Links.data(1,3)+1) + d.q;
end
    
for ii = 1:Nb.Links
    d.q = rho.cha(Links.data(ii,4)+1,Links.data(ii,5)+1,Links.data(ii,6)+1) *d.x*d.y*d.z;
    if (d.q>=0)
        Q.cha.plus(Links.data(ii,6)+1)  = Q.cha.plus(Links.data(ii,6)+1)  + d.q;
    elseif (d.q<=0)
        Q.cha.minus(Links.data(ii,6)+1) = Q.cha.minus(Links.data(ii,6)+1) + d.q;
    end
    
    alt_histogram.data(Links.data(ii,6)+1) = alt_histogram.data(Links.data(ii,6)+1) +1;
end
alt_histogram.max = max(alt_histogram.data);
z                 = ((0:N.z-1)*d.z+z_gnd);

%% Plots
figure(1);
plot(alt_histogram.data,z,'k-',-Q.cha.minus,z,'b',Q.cha.plus,z,'r');
axis([0 alt_histogram.max FocusArea.z(1) FocusArea.z(2)]);

box on
title('Altitude histogram',...
    'Units','normalized','FontName',Font.Name,'Fontsize',Font.Size.Labels)
text(alt_histogram.max,FocusArea.z(2),[int2str(Nb.Points),' Pts'],...
    'FontName',Font.Name,'Fontsize',Font.Size.Labels,...
    'VerticalAlignment','Top','HorizontalAlignment','Right')

xlabel('Q_{cha} (C),n ()','FontName',Font.Name,'Fontsize',Font.Size.Labels)
ylabel('z (km)','FontName',Font.Name,'Fontsize',Font.Size.Labels)
set(gca,'XMinorTick','on','YMinorTick','on')
grid on

format long
fprintf('Total charge induced         : %f C.\n',sum(sum(sum(rho.cha))));
fprintf('Total positive induced charge: %f C.\n',sum(Q.cha.plus));
fprintf('Total negative induced charge: %f C.\n',sum(Q.cha.minus));