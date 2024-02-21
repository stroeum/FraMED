function [q,CWC] = GZ_Electrification_Scheme(temp,P_dry,w,dim,params)
    %% Scaled temperature:
    temp_celcius = temp - 273.15;       % convert temperature to degrees celcius
    const_a = -1.7*(10^-5);
    const_b = -0.003;
    const_c = -0.05;
    const_d = 0.13;
    
    % Finding the scaled temperature at each gridpoint, Equation 17 in Mansell2005:
    tau_H2O = (-21/params.T_reversal_H2O).*(temp-params.T_triplet_H2O); % scaled temperature, unitless/degrees celcius
    f_tau_H2O = (const_a*(tau_H2O.^3))+(const_b*(tau_H2O.^2))+(const_c*tau_H2O)+const_d;
    tau_NH3 = (-21/params.T_reversal_NH3)*(temp-params.T_triplet_H2O); % scaled temperature, unitless/degrees celcius
    f_tau_NH3 = (const_a*(tau_NH3.^3))+(const_b*(tau_NH3.^2))+(const_c*tau_NH3)+const_d;
    
    % Finding delta L, pg. 6 of Mansell2005:
    CWC.crit_H2O = (temp.*0)+0.1; % critical cloud water content, g/m^3
    CWC.crit_NH3 = (temp.*0)+0.1; % critical cloud water content, g/m^3
    
    %% Calculate CWC in g/m^3 density:
    Ratmo = 3637; % gas constant for Jupiter's atmosphere, in J/(kg*K)
    rho_dry = P_dry./(Ratmo*temp); % density of the dry atmosphere (kg/m^3)
    
    CWC.ice_H2O = w.ice_H2O.*rho_dry*1000; % ice water content in g/m^3
    CWC.vapor_H2O = w.vapor_H2O.*rho_dry*1000; % vapor water content in g/m^3
    CWC.liquid_H2O = w.liquid_H2O.*rho_dry*1000; % liquid cloud water content in g/m^3
    CWC.rain_H2O = w.rain_H2O.*rho_dry*1000; % rain water content in g/m^3
    CWC.snow_H2O = w.snow_H2O.*rho_dry*1000; % snow water content in g/m^3
    CWC.all_H2O = CWC.ice_H2O + CWC.vapor_H2O + CWC.liquid_H2O + CWC.rain_H2O + CWC.snow_H2O; % CWC if we include everything
    CWC.only_liq_H2O = CWC.rain_H2O + CWC.liquid_H2O; % CWC if we include just water
    CWC.icesnow_H2O = CWC.ice_H2O + CWC.snow_H2O; % CWC if we include just ice and snow
    
    CWC.ice_NH3 = w.ice_NH3.*rho_dry*1000; % ice water content in g/m^3
    CWC.vapor_NH3 = w.vapor_NH3.*rho_dry*1000; % vapor water content in g/m^3
    CWC.liquid_NH3 = w.liquid_NH3.*rho_dry*1000; % liquid cloud water content in g/m^3
    CWC.rain_NH3 = w.rain_NH3.*rho_dry*1000; % rain water content in g/m^3
    CWC.snow_NH3 = w.snow_NH3.*rho_dry*1000; % snow water content in g/m^3
    CWC.all_NH3 = CWC.ice_NH3 + CWC.vapor_NH3 + CWC.liquid_NH3 + CWC.rain_NH3 + CWC.snow_NH3; % CWC if we include everything
    CWC.only_liq_NH3 = CWC.rain_NH3 + CWC.liquid_NH3; % CWC if we include just water
    CWC.icesnow_NH3 = CWC.ice_NH3 + CWC.snow_NH3; % CWC if we include just ice and snow
    
    totalWater = w.liquid_H2O+w.rain_H2O;
    [~,index3] = max(max(max(totalWater)));
    [~,index2] = max(max(totalWater(:,:,index3)));
    [~,index1] = max(totalWater(:,index2,index3));
    if max(max(max(abs(w.liquid_H2O+w.rain_H2O))))<(10^-6)
        q.total = 0.*w.rain_H2O;
        fprintf('\nNot enough liquid to apply electrification scheme.\n');
        return
    end
    % Finding n_ice and n_snow an easy way
    N_A = 6.022*(10^23); % Avogadro's number, particles/mol
    
    %% Finding delta Ls - different delta L for each species:
    dL_rain_H2O = temp.*0;
    dL_liquid_H2O = dL_rain_H2O;
    dL_rain_NH3 = dL_rain_H2O; 
    dL_liquid_NH3 = dL_rain_H2O; 

    % Based on Equation 16 in Mansell2005:
    for i = params.minLon:1:params.maxLon
        for j = params.minLat:1:params.maxLat
            for k = 1:1:dim.pressure
                if temp_celcius(i,j,k) > params.T_reversal_H2O && temp_celcius(i,j,k)<=0 
                    if max(w.rain_H2O(i,j,:)) >= (10^-6)
                        dL_rain_H2O(i,j,k) = CWC.rain_H2O(i,j,k)-CWC.crit_H2O(i,j,k);     % find dL in g/m^3
                    end
                    if max(w.liquid_H2O(i,j,:)) >= (10^-6)
                        dL_liquid_H2O(i,j,k) = CWC.liquid_H2O(i,j,k)-CWC.crit_H2O(i,j,k);     % find dL in g/m^3
                    end
                elseif temp_celcius(i,j,k) <= params.T_reversal_H2O
                    if max(w.rain_H2O(i,j,:)) >= (10^-6)
                        dL_rain_H2O(i,j,k) = CWC.rain_H2O(i,j,k);     % find dL in g/m^3
                    end
                    if max(w.liquid_H2O(i,j,:)) >= (10^-6)
                        dL_liquid_H2O(i,j,k) = CWC.liquid_H2O(i,j,k);     % find dL in g/m^3
                    end
                end
                if temp_celcius(i,j,k) > params.T_reversal_NH3 && temp_celcius(i,j,k)<=0 && (w.liquid_NH3(i,j,k)+w.rain_NH3(i,j,k)) >= (10^-6)
                    dL_rain_NH3(i,j,k) = CWC.rain_NH3(i,j,k)-CWC.crit_NH3(i,j,k);     % find dL in g/m^3
                    dL_liquid_NH3(i,j,k) = CWC.liquid_NH3(i,j,k)-CWC.crit_NH3(i,j,k);     % find dL in g/m^3
                elseif temp_celcius(i,j,k) <= params.T_reversal_NH3 && (w.liquid_NH3(i,j,k)+w.rain_NH3(i,j,k)) >= (10^-6)
                    dL_rain_NH3(i,j,k) = CWC.rain_NH3(i,j,k);     % find dL in g/m^3
                    dL_liquid_NH3(i,j,k) = CWC.liquid_NH3(i,j,k);     % find dL in g/m^3
                end
            end
        end
    end
    
    %% Table 1 in Palotai2008:
    % General values:
    P_0 = 100000;               % reference pressure, in Pa (1 bar)
    gamma_ice = 0.47;           % gamma is the same for snow and ice
    gamma_liquid = 0.33;        % gamma is the same for rain and liquid
    
    % For water:
    alpha_H2O = 7.06165*10^-3; % alpha constant for water ice
    beta_H2O = 2.0; % beta constant for water ice
    c_H2O = 5.38*10^7;
    d_H2O = 0.75;
    x_snow_H2O = 50.172;
    y_snow_H2O = 0.47798;
    x_rain_H2O = 2615.1;
    y_rain_H2O = 0.74245;
    rho_snow_H2O = 0.5*917;     % density of snow (half density of ice), kg/m^3;
    
    % For ammonia:
    alpha_NH3 = 7.06165*10^-3; % alpha constant for water ice
    beta_NH3 = 2.0; % beta constant for water ice
    c_NH3 = 5.38*10^7;
    d_NH3 = 0.75;
    x_snow_NH3 = 48.606;
    y_snow_NH3 = 0.488803;
    x_rain_NH3 = 2479;
    y_rain_NH3 = 0.76217;
    rho_snow_NH3 = 0.5*786.8;   % density of snow (half density of ice), kg/m^3;

    N_ice_H2O = temp.*0; % initialized params.number density of the ice, 1/m^3
    M_ice_H2O = N_ice_H2O; % initialized mass of the ice, kg
    D_ice_H2O = N_ice_H2O; % initialized diameter of the ice, m
    N_snow_H2O = temp.*0; % initialized params.number density of the ice, 1/m^3
    M_snow_H2O = N_snow_H2O; % initialized mass of the ice, kg
    D_snow_H2O = N_snow_H2O; % initialized diameter of the ice, m

    N_ice_NH3 = temp.*0; % initialized params.number density of the ice, 1/m^3
    M_ice_NH3 = N_ice_NH3; % initialized mass of the ice, kg
    D_ice_NH3 = N_ice_NH3; % initialized diameter of the ice, m
    N_snow_NH3 = temp.*0; % initialized params.number density of the ice, 1/m^3
    M_snow_NH3 = N_snow_NH3; % initialized mass of the ice, kg
    D_snow_NH3 = N_snow_NH3; % initialized diameter of the ice, m
    
    %% Initialize all mass weighted average terminal fall speeds
    v_rain_H2O = temp.*0;       % initialized mass-avg. term. velocity of rain, m/s
    v_liquid_H2O = temp.*0;     % initialized mass-avg. term. velocity of liquid, m/s
    v_snow_H2O = temp.*0;       % initialized mass-avg. term. velocity of snow, m/s
    v_ice_H2O = temp.*0;        % initialized mass-avg. term. velocity of ice, m/s
    N_0_snow_H2O = temp.*0;     % initialized intercept parameter
    
    v_rain_NH3 = temp.*0;       % initialized mass-avg. term. velocity of rain, m/s
    v_liquid_NH3 = temp.*0;     % initialized mass-avg. term. velocity of liquid, m/s
    v_snow_NH3 = temp.*0;       % initialized mass-avg. term. velocity of snow, m/s
    v_ice_NH3 = temp.*0;        % initialized mass-avg. term. velocity of ice, m/s
    N_0_snow_NH3 = temp.*0;     % initialized intercept parameter
    
    % Based on Equations 6-8 in Palotai2008
    N_ice_H2O = c_H2O*(rho_dry.* w.ice_H2O).^d_H2O;
    N_ice_NH3 = c_NH3*(rho_dry.* w.ice_NH3).^d_NH3;
    
    for i = params.minLon:1:params.maxLon
        for j = params.minLat:1:params.maxLat
            for k = 1:1:dim.pressure
                if N_ice_H2O(i,j,k) ~= 0 % avoid dividing by zero
                    M_ice_H2O(i,j,k) = (rho_dry(i,j,k)*w.ice_H2O(i,j,k))/N_ice_H2O(i,j,k);
                    D_ice_H2O(i,j,k) = (M_ice_H2O(i,j,k)/alpha_H2O)^(1/beta_H2O);
                    v_ice_H2O(i,j,k) = x_snow_H2O*(D_ice_H2O(i,j,k)^y_snow_H2O);
                end
                if N_ice_NH3(i,j,k) ~= 0 % avoid dividing by zero
                    M_ice_NH3(i,j,k) = (rho_dry(i,j,k)*w.ice_NH3(i,j,k))/N_ice_NH3(i,j,k);
                    D_ice_NH3(i,j,k) = (M_ice_NH3(i,j,k)/alpha_NH3)^(1/beta_NH3);
                    v_ice_NH3(i,j,k) = x_snow_NH3*(D_ice_NH3(i,j,k)^y_snow_NH3);
                end
            end
        end
    end
    
    % Determining N_0_rain;
    N_0_rain_H2O = 8*(10^6); % intercept parameter for rain, pg. 24 Palotai2008
    rho_rain_H2O = 1000.0; %density of cloud water and rain, kg/m^3
    lambda_rain_H2O = temp.*0;
    lambda_snow_H2O = lambda_rain_H2O;
    lambda_liquid_H2O = lambda_rain_H2O;
    N_0_rain_NH3 = 8*(10^6); % intercept parameter for rain, pg. 24 Palotai2008
    rho_rain_water_NH3 = 733.0; %density of cloud water and rain, kg/m^3
    lambda_rain_NH3 = temp.*0;
    lambda_snow_NH3 = lambda_rain_NH3;
    lambda_liquid_NH3 = lambda_rain_NH3;
    
    D_min = 200*(10^-6); % minimum raindrop diameter, meters
    D_max = 5*(10^-3); % maximum raindrop diameter, meters
    rain_diameters = linspace(D_min, D_max, 100);
    snow_diameters = linspace((500*10^-6),(5*10^-3),100); % snow diameters, pg. 7 of Palotai2008
    liquid_diameters = linspace((10^-7),(200*10^-6),100); % cloud water droplet diameters
    ice_diameters = linspace((10^-6),(500*10^-6),100); % ice diameters, pg. 8 of Palotai2008
    D_snow_H2O = mean(snow_diameters);
    D_snow_NH3 = D_snow_H2O;
    
    N_rain_H2O = temp.*0; % particle density, 1/m^3
    N_snow_H2O = N_rain_H2O; 
    N_liquid_H2O = N_rain_H2O; 
    N_rain_NH3 = temp.*0; % particle density, 1/m^3
    N_snow_NH3 = N_rain_NH3; 
    N_liquid_NH3 = N_rain_NH3; 
    
    % Based on Equations 9-14 and 28 in Palotai2008
    for i = params.minLon:1:params.maxLon
        for j = params.minLat:1:params.maxLat
            for k = 1:1:dim.pressure
            % H2O Related:
                % RAIN:
                if w.rain_H2O(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_rain_H2O(i,j,k) = ((pi*rho_rain_H2O*N_0_rain_H2O)/(rho_dry(i,j,k)*w.rain_H2O(i,j,k)))^(1/4);
                    N_rain_H2O(i,j,k) = trapz(rain_diameters,N_0_rain_H2O*exp(-lambda_rain_H2O(i,j,k)*rain_diameters));
                    v_rain_H2O(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_liquid)*(x_rain_H2O/(lambda_rain_H2O(i,j,k)^y_rain_H2O))*(gamma(y_rain_H2O+4)/6);
                end
                % LIQUID: 
                if w.liquid_H2O(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_liquid_H2O(i,j,k) = ((pi*rho_rain_H2O*N_0_rain_H2O)/(rho_dry(i,j,k)*w.liquid_H2O(i,j,k)))^(1/4);
                    N_liquid_H2O(i,j,k) = trapz(liquid_diameters,N_0_rain_H2O*exp(-lambda_liquid_H2O(i,j,k)*liquid_diameters));
                    v_liquid_H2O(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_liquid)*(x_rain_H2O/(lambda_liquid_H2O(i,j,k)^y_rain_H2O))*(gamma(y_rain_H2O+4)/6);
                end
                % SNOW:
                if 0.01*exp(-0.12*(temp(i,j,k)-params.T_triplet_H2O))<1
                    N_0_snow_H2O(i,j,k) = 2*(10^8)*0.01*exp(-0.12*(temp(i,j,k)-params.T_triplet_H2O));
                else
                    N_0_snow_H2O(i,j,k) = 2*(10^8);
                end
                if w.snow_H2O(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_snow_H2O(i,j,k) = ((pi*rho_snow_H2O*N_0_snow_H2O(i,j,k))/(rho_dry(i,j,k)*w.snow_H2O(i,j,k)))^(1/4);
                    N_snow_H2O(i,j,k) = trapz(snow_diameters,N_0_snow_H2O(i,j,k)*exp(-lambda_snow_H2O(i,j,k)*snow_diameters));
                    v_snow_H2O(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_ice)*(x_snow_H2O/(lambda_snow_H2O(i,j,k)^y_snow_H2O))*(gamma(y_snow_H2O+4)/6);
                end
                % ICE:

    
            % NH3 Related:
                % RAIN:
                if w.rain_NH3(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_rain_NH3(i,j,k) = ((pi*rho_rain_water_NH3*N_0_rain_NH3)/(rho_dry(i,j,k)*w.rain_NH3(i,j,k)))^(1/4);
                    N_rain_NH3(i,j,k) = trapz(rain_diameters,N_0_rain_NH3*exp(-lambda_rain_NH3(i,j,k)*rain_diameters));
                    v_rain_NH3(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_liquid)*(x_rain_NH3/(lambda_rain_NH3(i,j,k)^y_rain_NH3))*(gamma(y_rain_NH3+4)/6);
                end
                % LIQUID:
                if w.liquid_NH3(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_liquid_NH3(i,j,k) = ((pi*rho_rain_water_NH3*N_0_rain_NH3)/(rho_dry(i,j,k)*w.liquid_NH3(i,j,k)))^(1/4);
                    N_liquid_NH3(i,j,k) = trapz(liquid_diameters,N_0_rain_NH3*exp(-lambda_liquid_NH3(i,j,k)*liquid_diameters));
                    v_liquid_NH3(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_liquid)*(x_rain_NH3/(lambda_liquid_NH3(i,j,k)^y_rain_NH3))*(gamma(y_rain_NH3+4)/6);
                end
                % SNOW:
                if 0.01*exp(-0.12*(temp(i,j,k)-params.T_triplet_NH3))<1
                    N_0_snow_NH3(i,j,k) = 2*(10^8)*0.01*exp(-0.12*(temp(i,j,k)-params.T_triplet_NH3));
                else
                    N_0_snow_NH3(i,j,k) = 2*(10^8);
                end
                if w.snow_NH3(i,j,k)*rho_dry(i,j,k) ~= 0 % avoid dividing by zero
                    lambda_snow_NH3(i,j,k) = ((pi*rho_snow_NH3*N_0_snow_NH3(i,j,k))/(rho_dry(i,j,k)*w.snow_NH3(i,j,k)))^(1/4);
                    N_snow_NH3(i,j,k) = trapz(snow_diameters,N_0_snow_NH3(i,j,k)*exp(-lambda_snow_NH3(i,j,k)*snow_diameters));
                    v_snow_NH3(i,j,k) = ((P_0/P_dry(i,j,k))^gamma_ice)*(x_snow_NH3/(lambda_snow_NH3(i,j,k)^y_snow_NH3))*(gamma(y_snow_NH3+4)/6);
                end
            end
        end
    end
    
    %% Determine charge separation of each species:
    % Based on Equation 14 in Mansell2005 / Equation 5 in Gardiner1985:
    dq_rain_ice_H2O = 7.3*((D_ice_H2O).^4).*((v_ice_H2O-v_rain_H2O).^3).*dL_rain_H2O.*f_tau_H2O;         % collisions of rain & ice
    dq_rain_snow_H2O = 7.3*((D_snow_H2O).^4).*((v_snow_H2O-v_rain_H2O).^3).*dL_rain_H2O.*f_tau_H2O;      % collisions of rain & snow

    dq_liquid_ice_H2O = 7.3*((D_ice_H2O).^4).*((v_ice_H2O-v_liquid_H2O).^3).*dL_liquid_H2O.*f_tau_H2O;     % collisions of liquid & ice
    dq_liquid_snow_H2O = 7.3*((D_snow_H2O).^4).*((v_snow_H2O-v_liquid_H2O).^3).*dL_liquid_H2O.*f_tau_H2O;  % collisions of liquid & snow
    
    dq_rain_ice_NH3 = 7.3*((D_ice_NH3).^4).*((v_ice_NH3-v_rain_NH3).^3).*dL_rain_NH3.*f_tau_NH3;         % collisions of rain & ice
    dq_rain_snow_NH3 = 7.3*((D_snow_NH3).^4).*((v_snow_NH3-v_rain_NH3).^3).*dL_rain_NH3.*f_tau_NH3;      % collisions of rain & snow

    dq_liquid_ice_NH3 = 7.3*((D_ice_NH3).^4).*((v_ice_NH3-v_liquid_NH3).^3).*dL_liquid_NH3.*f_tau_NH3;     % collisions of liquid & ice
    dq_liquid_snow_NH3 = 7.3*((D_snow_NH3).^4).*((v_snow_NH3-v_liquid_NH3).^3).*dL_liquid_NH3.*f_tau_NH3;  % collisions of liquid & snow
    
    % Assign signs to the charge separations:
    % Based on Section 5.4 in Ziegler1991:
    for i = params.minLon:1:params.maxLon
        for j = params.minLat:1:params.maxLat
            for k = 1:1:dim.pressure
                if CWC.only_liq_H2O(i,j,k) > CWC.crit_H2O(i,j,k) && temp_celcius(i,j,k)>params.T_reversal_H2O
                    dq_rain_ice_H2O(i,j,k) = abs(dq_rain_ice_H2O(i,j,k));
                    dq_rain_snow_H2O(i,j,k) = abs(dq_rain_snow_H2O(i,j,k));
                    dq_liquid_ice_H2O(i,j,k) = abs(dq_liquid_ice_H2O(i,j,k));
                    dq_liquid_snow_H2O(i,j,k) = abs(dq_liquid_snow_H2O(i,j,k));
                else
                    dq_rain_ice_H2O(i,j,k) = -1*abs(dq_rain_ice_H2O(i,j,k));
                    dq_rain_snow_H2O(i,j,k) = -1*abs(dq_rain_snow_H2O(i,j,k));
                    dq_liquid_ice_H2O(i,j,k) = -1*abs(dq_liquid_ice_H2O(i,j,k));
                    dq_liquid_snow_H2O(i,j,k) = -1*abs(dq_liquid_snow_H2O(i,j,k));
                end
                
                if CWC.only_liq_NH3(i,j,k) > CWC.crit_NH3(i,j,k) && temp_celcius(i,j,k)>params.T_reversal_NH3
                    dq_rain_ice_NH3(i,j,k) = abs(dq_rain_ice_NH3(i,j,k));
                    dq_rain_snow_NH3(i,j,k) = abs(dq_rain_snow_NH3(i,j,k));
                    dq_liquid_ice_NH3(i,j,k) = abs(dq_liquid_ice_NH3(i,j,k));
                    dq_liquid_snow_NH3(i,j,k) = abs(dq_liquid_snow_NH3(i,j,k));
                else
                    dq_rain_ice_NH3(i,j,k) = -1*abs(dq_rain_ice_NH3(i,j,k));
                    dq_rain_snow_NH3(i,j,k) = -1*abs(dq_rain_snow_NH3(i,j,k));
                    dq_liquid_ice_NH3(i,j,k) = -1*abs(dq_liquid_ice_NH3(i,j,k));
                    dq_liquid_snow_NH3(i,j,k) = -1*abs(dq_liquid_snow_NH3(i,j,k));
                end
            end
        end
    end
    
    % Find the total charge density:
    q.rain_ice_H2O = (N_ice_H2O.*dq_rain_ice_H2O);
    q.rain_snow_H2O = (N_snow_H2O.*dq_rain_snow_H2O);
    q.liquid_ice_H2O = (N_ice_H2O.*dq_liquid_ice_H2O);
    q.liquid_snow_H2O = (N_snow_H2O.*dq_liquid_snow_H2O);
    q.total_H2O = q.rain_ice_H2O + q.rain_snow_H2O + q.liquid_ice_H2O + q.liquid_snow_H2O;

    q.rain_ice_NH3 = (N_ice_NH3.*dq_rain_ice_NH3);
    q.rain_snow_NH3 = (N_snow_NH3.*dq_rain_snow_NH3);
    q.liquid_ice_NH3 = (N_ice_NH3.*dq_liquid_ice_NH3);
    q.liquid_snow_NH3 = (N_snow_NH3.*dq_liquid_snow_NH3);
    q.total_NH3 = q.rain_ice_NH3 + q.rain_snow_NH3 + q.liquid_ice_NH3 + q.liquid_snow_NH3;

    q.total = q.total_H2O + q.total_NH3;
    fprintf('\nGZ Scheme script successfully executed.\n');
end