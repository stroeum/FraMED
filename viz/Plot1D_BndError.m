close all
clearvars -except sims

beep  off
cd ../results
load BndError.dat
cd ../viz

A1 = BndError(:,1);
A2 = BndError(:,2);
A3 = BndError(:,3);

A = abs((A1-A2)./(A2+A3));
plot(A*100)
xlabel('steps','FontSize',12)
ylabel('\epsilon_{BC}','FontSize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
grid on