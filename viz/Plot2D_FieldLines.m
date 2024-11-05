function [] = Plot2D_FieldLines(num,sims)
    fprintf('\n*** Executing Plot2D_FieldLines.m function. ***\n');
    if ~exist('sims','var') || ~isfield(sims,'pathPNGs') 
        prompt1 = "\nWhat is the planetary body that the simulation is focused on? (No quotation marks needed for string input)\n-->";
        sims.objectName = input(prompt1,'s');
        prompt2 = "\nWhat type of discharge is this? (Leader / Streamer)\n-->";
        sims.objectType = input(prompt2,'s');
        while ~strcmp(sims.objectType,'Streamer') && ~strcmp(sims.objectType,'Leader')
            fprintf('\n\tNot an acceptable input. Please enter Streamer or Leader.\n');
            sims.objectType = input(prompt2,'s');
        end
    
        % Settings to ensure proper directory referencing:
        sims.pathPNGs = ['../Figures/',sims.objectName,'/',sims.objectType,'/PNGs'];
        if ~exist(sims.pathPNGs,'dir')
            mkdir(sims.pathPNGs);
        end
    end 

    cd ../results
    dxyz         = load('dxyz.dat');
    Nxyz         = load('Nxyz.dat');
    
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
    
    z_gnd        = load('z_gnd.dat');
    ChargeLayers = load('ChargeLayers.dat');
    cd ../viz
    dx = dxyz(1);               % _m
    dy = dxyz(2);               % _m
    dz = dxyz(3);               % _m
    
    Nx = Nxyz(1);
    Ny = Nxyz(2);
    Nz = Nxyz(3);
    
    Lx = (Nx-1)*dx*1e-3;         % _m
    Ly = (Ny-1)*dy*1e-3;         % _m
    Lz = (Nz-1)*dz*1e-3;         % _m
    
    clear dxyz Nxyz
    
    fprintf('\tData loaded\n')
    
    %% Fix plot parameters
    NbChargeLayers = size(ChargeLayers);
    NbChargeLayers = NbChargeLayers(1);
    ChargeLayersLineStyle  = '-';
    ChargeLayersLineWidth  = 1;
    
    %% Derive data for plotting
    x        = (0:Nx-1)'*dx*1e-3;
    y        = (0:Ny-1)'*dy*1e-3;
    z        = ((0:Nz-1)'*dz+z_gnd)*1e-3;
    [y,z1]    = meshgrid(y,z);
    [x,z2]    = meshgrid(x,z);
    
    %% Plot
    figure;
    hold on
    %set(gcf,'Units','inches','OuterPosition', [10 10 20 30]/3)
    axis([min(y(:)) max(y(:)) min(z1(:)) max(z1(:))]);
    streamslice(y,z1,Ey2D,Ez2D,10,'arrows');
    % quiver(y,z,Ey2D',Ez2D');
    set(findobj('Type','line'),'Color','k')
    plot([min(y(:)) max(y(:))], [z_gnd z_gnd]*1e-3,'k');
    for ii=1:NbChargeLayers
        if ChargeLayers(ii,1)>0
            tempColor = 'r'; % positive charge region
        else
            tempColor = 'b'; % negative charge region
        end
        rectangle('Position',[(ChargeLayers(ii,3)-2*ChargeLayers(ii,6)/2)*1e-3,(z_gnd+ChargeLayers(ii,4)-ChargeLayers(ii,7)/2)*1e-3,2*ChargeLayers(ii,6)*1e-3,ChargeLayers(ii,7)*1e-3],...
            'Curvature',[0,0],...
            'LineWidth',ChargeLayersLineWidth,'LineStyle',ChargeLayersLineStyle,'EdgeColor',tempColor);
    %     text((ChargeLayers(ii,3)+ChargeLayers(ii,6))*1e-3,(z_gnd+ChargeLayers(ii,4))*1e-3,...
    %         ['\leftarrow',num2str(ChargeLayers(ii,1),3),' C'],...
    %         'HorizontalAlignment','left','BackgroundColor','w',...
    %         'FontSize',10,'Color',Color(ii))
    end
    
    hold off
    xlabel('$y$ (km)','Interpreter','latex','FontSize',16);
    ylabel('$z$ (km)','Interpreter','latex','FontSize',16);
    set(gca,'FontSize',10,'TickLabelInterpreter','latex');
    %title(['Field Lines on ',sims.objectName],'Interpreter','latex','FontSize',18);
        
    axis xy
    axis image
    box on
    
    exportgraphics(gcf,[sims.pathPNGs,'/FieldLines_',sims.objectName,'_',sims.objectType,'.png'],'BackgroundColor','white','Resolution',300);
end 