function exportChargeDensity(params,Nxyz,FraMED)
    fileID = fopen([params.path2datafiles,params.discharge_type,'_',params.typeWC,'_rhoAmb.dat'],'w');
    for kk = 1:Nxyz(3)
        for ii = 1:Nxyz(1)
            for jj = 1:Nxyz(2)
                fprintf(fileID,'%.6f ',(params.SourcesCPP_factor*FraMED.rho(kk)));
            end
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end