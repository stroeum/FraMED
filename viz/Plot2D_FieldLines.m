function [] = Plot2D_FieldLines(num,sims)
    fprintf('\n*** Executing Plot2D_FieldLines.m function. ***\n');
    if ~exist('sims','var') 
        specifySimDetails;
    end 

    cd ../results
    
    if ~exist('num','var') || num=='~'
        Ex2D         = load('Ex2d.dat');
        Ey2D         = load('Ey2d.dat');
        Ez2D         = load('Ez2d.dat');
        phi2D        = load('phi2d.dat');
    else 
        Ex2D         = load(['Ex2D',num2str(num),'.dat']);
        Ey2D         = load(['Ey2D',num2str(num),'.dat']);
        Ez2D         = load(['Ez2D',num2str(num),'.dat']);
        phi2D        = load(['phi2D',num2str(num),'.dat']);
    end
    
    ChargeLayers = load('ChargeLayers.dat');
    cd ../viz
    
    fprintf('\tData loaded\n')
    
    %% Fix plot parameters
    NbChargeLayers = size(ChargeLayers);
    NbChargeLayers = NbChargeLayers(1);
    ChargeLayersLineStyle  = '-';
    ChargeLayersLineWidth  = 1;
    
    %% Derive data for plotting
    x        = (0:sims.domain.dx:sims.domain.maxx)';
    y        = (0:sims.domain.dy:sims.domain.maxy)';
    z        = (sims.domain.gnd:sims.domain.dz:sims.domain.maxz)';
    [y,z1]    = meshgrid(y,z);
    [x,z2]    = meshgrid(x,z);
    
    %% Plot
    figure;
    hold on
    %set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/3)
    axis([sims.spatialFactor.Number*min(y(:)) sims.spatialFactor.Number*max(y(:)) sims.spatialFactor.Number*min(z1(:)) sims.spatialFactor.Number*max(z1(:))]);
    streamslice(sims.spatialFactor.Number*y,sims.spatialFactor.Number*z1,Ey2D,Ez2D,10,'arrows');
    % quiver(y,z,Ey2D',Ez2D');
    set(findobj('Type','line'),'Color','k')
    plot([sims.spatialFactor.Number*min(y(:)) sims.spatialFactor.Number*max(y(:))], [sims.domain.gnd sims.domain.gnd]*sims.spatialFactor.Number,'k');
    for ii=1:NbChargeLayers
        if ChargeLayers(ii,1)>0
            tempColor = 'r'; % positive charge region
        else
            tempColor = 'b'; % negative charge region
        end
        rectangle('Position',[(ChargeLayers(ii,3)-2*ChargeLayers(ii,6)/2)*sims.spatialFactor.Number,(sims.domain.gnd+ChargeLayers(ii,4)-ChargeLayers(ii,7)/2)*sims.spatialFactor.Number,2*ChargeLayers(ii,6)*sims.spatialFactor.Number,ChargeLayers(ii,7)*sims.spatialFactor.Number],...
            'Curvature',[0,0],...
            'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',tempColor);
    %     text((ChargeLayers(ii,3)+ChargeLayers(ii,6))*sims.spatialFactor.Number,(sims.domain.gnd+ChargeLayers(ii,4))*sims.spatialFactor.Number,...
    %         ['\leftarrow',num2str(ChargeLayers(ii,1),3),' C'],...
    %         'HorizontalAlignment','left','BackgroundColor','w',...
    %         'FontSize',10,'Color',Color(ii))
    end
    
    hold off
    xlabel(strcat('$y$ (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',16);
    ylabel(strcat('$z$ (',sims.spatialFactor.LaTeX,'m)'),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',10,'TickLabelInterpreter','latex');
    %title(strcat('Field Lines on ',sims.objectName),'Interpreter','latex','FontSize',18);
        
    axis xy
    axis image
    box on
    
    exportgraphics(gcf,strcat(sims.pathPNGs,'/FieldLines_',sims.objectName,'_',sims.objectType,'.png'),'BackgroundColor','white','Resolution',300);
end 