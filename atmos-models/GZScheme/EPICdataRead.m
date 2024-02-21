function [temp,lat,lon,w,P_dry,P] = EPICdataRead(params)
    % Read in data from EPIC output .nc:
    temp = ncread([params.basefolder,params.filename],'t');          % temperature, Kelvin
    lat = ncread([params.basefolder,params.filename],'lat_u');       % latitude, degrees
    lon = ncread([params.basefolder,params.filename],'lon_u');       % longitude, degrees
    if length(lon)>1 || length(lat)>1
        fprintf(['The file below has multiple latitude/longitude values.\n\tLatitudes: ',num2str(length(lat)),'\n\tLongitudes: ',num2str(length(lon))]);
    end
    
    % Cloud water mixing ratio:
    w.ice_H2O = ncread([params.basefolder,params.filename],'H_2O_solid');     % water ice mixing ratio, kg/kg
    w.vapor_H2O = ncread([params.basefolder,params.filename],'H_2O_vapor');   % water vapor mixing ratio, kg/kg
    w.liquid_H2O = ncread([params.basefolder,params.filename],'H_2O_liquid'); % water liquid cloud mixing ratio, kg/kg
    w.rain_H2O = ncread([params.basefolder,params.filename],'H_2O_rain');     % water rain mixing ratio, kg/kg
    w.snow_H2O = ncread([params.basefolder,params.filename],'H_2O_snow');     % snow mixing ratio, kg/kg
    
    % Cloud ammonia mixing ratio:
    w.ice_NH3 = ncread([params.basefolder,params.filename],'NH_3_solid');     % water ice mixing ratio, kg/kg
    w.vapor_NH3 = ncread([params.basefolder,params.filename],'NH_3_vapor');   % water vapor mixing ratio, kg/kg
    w.liquid_NH3 = ncread([params.basefolder,params.filename],'NH_3_liquid'); % water liquid cloud mixing ratio, kg/kg
    w.rain_NH3 = ncread([params.basefolder,params.filename],'NH_3_rain');     % water rain mixing ratio, kg/kg
    w.snow_NH3 = ncread([params.basefolder,params.filename],'NH_3_snow');     % snow mixing ratio, kg/kg
    
    P_dry = ncread([params.basefolder,params.filename],'pdry');           % pressure of dry atmoshpere, in Pa
    P = ncread([params.basefolder,params.filename],'p');                  % pressure, Pa
end