% LMA charges %
if      strcmp(Layers.Type,'disks')
    for ii=1:Nb.Layers
        subplot(h.xy.fig);
        rectangle('Position',[(Layers.data(ii,2)-Layers.data(ii,5)),(Layers.data(ii,3)-Layers.data(ii,6)),2*Layers.data(ii,5),2*Layers.data(ii,6)],'Curvature',[1,1],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
        subplot(h.yz.fig);
        rectangle('Position',[(z_gnd+Layers.data(ii,4)-Layers.data(ii,7)/2),(Layers.data(ii,3)-Layers.data(ii,6)),Layers.data(ii,7),2*Layers.data(ii,6)],'Curvature',[0,0],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
        subplot(h.xz.fig);
        rectangle('Position',[(Layers.data(ii,2)-Layers.data(ii,5)),(z_gnd+Layers.data(ii,4)-Layers.data(ii,7)/2),2*Layers.data(ii,5),Layers.data(ii,7)],'Curvature',[0,0],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
    end
elseif  strcmp(Layers.Type,'spheres')
    for ii=1:Nb.Layers
        subplot(h.xy.fig);
        rectangle('Position',[(Layers.data(ii,2)-Layers.data(ii,5)),(Layers.data(ii,3)-Layers.data(ii,6)),2*Layers.data(ii,5),2*Layers.data(ii,6)],'Curvature',[1,1],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
        subplot(h.yz.fig);
        rectangle('Position',[(z_gnd+Layers.data(ii,4)-Layers.data(ii,7)),(Layers.data(ii,3)-Layers.data(ii,6)),2*Layers.data(ii,7),2*Layers.data(ii,6)],'Curvature',[1,1],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
        subplot(h.xz.fig);
        rectangle('Position',[(Layers.data(ii,2)-Layers.data(ii,5)),(z_gnd+Layers.data(ii,4)-Layers.data(ii,7)),2*Layers.data(ii,5),2*Layers.data(ii,7)],'Curvature',[1,1],'LineWidth',Layers.Line.Width,'LineStyle',Layers.Line.Style,'EdgeColor',Layers.Edge.Color(ii,:));
    end
elseif  strcmp(Layers.Type,'arbitrary')
    subplot(h.yz.fig);
    hold on;
    contourf(z,y,rho.YZ,60,'LineColor','none');
    caxis([-max(max(abs(rho.YZ))) max(max(abs(rho.YZ)))])
    % imagesc(Y,Z,rho.YZ');
    colorbar('Location','East');
    hold off
    
    subplot(h.xz.fig);
    hold on
    contourf(x,z,rho.XZ',60,'LineColor','none');
    caxis([-max(max(abs(rho.XZ))) max(max(abs(rho.XZ)))])
    % imagesc(Y,Z,rho.XZ');
    colorbar('Location','North');
    hold off
end