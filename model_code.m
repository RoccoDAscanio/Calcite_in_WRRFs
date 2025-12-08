%% ------------- WRRF Model Code --------------

clc;
clear;

% ---------------- RUN FLAG -------------------

% flag options: 
% 'time_integrate' -> Solve for the time-integrated solution
% 'steady_state'   -> Solve for the steady state soluiton
run_flag = 'time_integrate';                                                % flag to determine how model is run

% ---------------- INPUT PARAMAETERS ----------

Flow_rate_1         = 33;                                                   % Plant flow (MGD)
Flow_rate           = Flow_rate_1*1E6*0.00378;                              % Plant flow (m3 / day) 
frac_WAS_flow       = 0.0075;                                               % Fraction of total flow removed in WAS (dimensionless)
ratio_WAS_solids    = 3;                                                    % Ratio of WAS TSS : Biological Tanks TSS (dimensionless)
CaCO3_in            = 200;                                                  % Concentration of CaCO3 dosed continulously to process in inflow (mg CaCO3 / L)
BOD_in              = 0.3;                                                  % BOD in influent (g / L)
BOD_COD             = 0.8;                                                  % BOD:COD ratio (dimensionless)
DIC_in              = 4000;                                                 % DIC concentration in influent (μmol / L)
ALK_in              = 3000;                                                 % Alkalinity in influent (μmol / L)
Ca_in               = 1000;                                                 % Calcium concentration in influent (umol / L)
R_carbon_1          = BOD_in.*(1/3.33)*(44.01/12).*1E3.*(1/44).*(1/BOD_COD);% CO2 generation rate (mol / m3)
R_carbon            = Flow_rate.*R_carbon_1;                                % CO2 generation rate (mol / day)
Depth_average       = 5;                                                    % Average depth of wastewater tanks (m)
Time_residence      = 1;                                                    % Residence time (day) 
T   = 25;                                                                   % Water temp (deg C) 
P   = 1;                                                                    % Pressure (atm)
S   = 1;                                                                    % Salinity (ppt)

n_power_diss     = 1;                                                       % order of reaction for dissolution (dimensionless) 
n_power_prec     = 1.7;                                                     % order of reaction for precipitation (dimensionless) 

k_mineral           = 1E-6;                                                 % Calcite dissolution rate constate (mol / m2 /s) 
n_all            = 0.0438;                                                  % Rction order WRT [H+] (dimensionless)
Calcium_activity = 1;                                                       % Fixed activity coefficient for Ca (dimensionless) 

% ---------- ODE parameters -------------- 

tmax     = 100;                                                             % final time for integration (day)
maxstep  = 0.1;                                                             % maximum step for integration (day)

% ---------- Calibration Parameters ------ 

alpha_ebull = 8.56300000000000;                                             % Ebullition factor (dimensionless) 
Dni_half_sat =  197.5000;                                                   % Denitrification half saturation (μM)
k_nitrif = 28.1712;                                                         % Nitirification Rate Constant (1/μM/yr)

% -------------------------------------------------------------------------
% ---------------- CaCO3 size distribution and rate constant --------------
% -------------------------------------------------------------------------

mineral_mass        = 100;                                                  % mass for caluclating size distribition (grams)
r_grain             = 'PSD';                                                

E_a                 = 13.72E3;                                              % apparent activation energy (J/mol)
R_gas               = 8.314;                                                % gas constant (J/K/mol)
T_K                = T + 273;                                               % Temperature (kelvin)

r_mineral = logspace(log10(1), log10(1000), 14);                            % 40 bins

% Gaussian parameters (in log10 space)
mu = log10(50);                                                             % peak at 70 μm
sigma = 0.25;                                                               % adjust for desired spread (narrower = smaller sigma)

% Compute Gaussian in log-space
log_r = log10(r_mineral);
K_mineral = exp(-((log_r - mu).^2) / (2 * sigma^2));

% Normalize to percent
K_mineral_percent = K_mineral / sum(K_mineral) * 100;

if strcmp(r_grain,'PSD')
    % ---------------------------------------------------------------------
    % --- fraction of grains in each radius bin [dimensionless] -----------
    % ---------------------------------------------------------------------
    % NOTE: distribution here follows that of Renforth [2012] | doi:10.1016/j.ijggc.2012.06.011

    K_mineral_1 =   K_mineral(1);   
    K_mineral_2 =   K_mineral(2);   
    K_mineral_3 =   K_mineral(3);   
    K_mineral_4 =   K_mineral(4);   
    K_mineral_5 =   K_mineral(5);   
    K_mineral_6 =   K_mineral(6);   
    K_mineral_7 =   K_mineral(7);   
    K_mineral_8 =   K_mineral(8);   
    K_mineral_9 =   K_mineral(9);   
    K_mineral_10 =  K_mineral(10);   
    K_mineral_11 =  K_mineral(11);  
    K_mineral_12 =  K_mineral(12);   
    K_mineral_13 =  K_mineral(13);   
    K_mineral_14 =  K_mineral(14);

    % ---------------------------------------------------------------------
    % --- grain radii [um] ------------------------------------------------
    % ---------------------------------------------------------------------
    r_mineral_1  =  r_mineral(1);
    r_mineral_2  =  r_mineral(2);
    r_mineral_3  =  r_mineral(3);
    r_mineral_4  =  r_mineral(4);
    r_mineral_5  =  r_mineral(5);
    r_mineral_6  =  r_mineral(6);
    r_mineral_7  =  r_mineral(7);
    r_mineral_8  =  r_mineral(8);
    r_mineral_9  =  r_mineral(9);    
    r_mineral_10 =  r_mineral(10);
    r_mineral_11 =  r_mineral(11);
    r_mineral_12 =  r_mineral(12);
    r_mineral_13 =  r_mineral(13);
    r_mineral_14 =  r_mineral(14);
else
    K_mineral = 1.0;
    r_mineral = r_grain;
end

% -------------------------------------------------------------------------
% --- calculate SSA and dissolution rate constant values ------------------
% -------------------------------------------------------------------------
% NOTE: SSA in units of m2/g
% NOTE: SSA calculated based on grain size according to original arrays
if strcmp(r_grain,'PSD')

    SSA_1             = (15.53*(r_mineral_1^-0.2829))  - 1.437;
    SSA_2             = (15.53*(r_mineral_2^-0.2829))  - 1.437;
    SSA_3             = (15.53*(r_mineral_3^-0.2829))  - 1.437;
    SSA_4             = (15.53*(r_mineral_4^-0.2829))  - 1.437;
    SSA_5             = (15.53*(r_mineral_5^-0.2829))  - 1.437;
    SSA_6             = (15.53*(r_mineral_6^-0.2829))  - 1.437;
    SSA_7             = (15.53*(r_mineral_7^-0.2829))  - 1.437;
    SSA_8             = (15.53*(r_mineral_8^-0.2829))  - 1.437;
    SSA_9             = (15.53*(r_mineral_9^-0.2829))  - 1.437;
    SSA_10            = (15.53*(r_mineral_10^-0.2829)) - 1.437;
    SSA_11            = (15.53*(r_mineral_11^-0.2829)) - 1.437;
    SSA_12            = (15.53*(r_mineral_12^-0.2829)) - 1.437;
    SSA_13            = (15.53*(r_mineral_13^-0.2829)) - 1.437;
    SSA_14            = (15.53*(r_mineral_14^-0.2829)) - 1.437;
else
    SSA               = (15.53*(r_mineral^-0.2829)) - 1.437;
end
% -------------------------------------------------------------------------
% NOTE: rate constants in units of 1/day
if strcmp(r_grain,'PSD')
    k_diss_mineral_1  = k_mineral.*SSA_1.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_2  = k_mineral.*SSA_2.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_3  = k_mineral.*SSA_3.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_4  = k_mineral.*SSA_4.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_5  = k_mineral.*SSA_5.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_6  = k_mineral.*SSA_6.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_7  = k_mineral.*SSA_7.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_8  = k_mineral.*SSA_8.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_9  = k_mineral.*SSA_9.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_10 = k_mineral.*SSA_10.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_11 = k_mineral.*SSA_11.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_12 = k_mineral.*SSA_12.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_13 = k_mineral.*SSA_13.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;
    k_diss_mineral_14 = k_mineral.*SSA_14.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;

    k_calcite_dis         = K_mineral_1.*k_diss_mineral_1 + K_mineral_2.*k_diss_mineral_2 + K_mineral_3.*k_diss_mineral_3 + ...
                      K_mineral_4.*k_diss_mineral_4 + K_mineral_5.*k_diss_mineral_5 + K_mineral_6.*k_diss_mineral_6 + ...
                      K_mineral_7.*k_diss_mineral_7 + K_mineral_8.*k_diss_mineral_8 + K_mineral_9.*k_diss_mineral_9 + ...
                      K_mineral_10.*k_diss_mineral_10 + K_mineral_11.*k_diss_mineral_11 + K_mineral_12.*k_diss_mineral_12 + ...
                      K_mineral_13.*k_diss_mineral_13 + K_mineral_14.*k_diss_mineral_14;  %1/day

else
        k_diss_mineral    = k_mineral.*SSA.*mineral_mass.*exp(-(E_a)./(R_gas.*T_K)).*60.*60.*24.*1E4;

        k_calcite_dis         = K_mineral.*k_diss_mineral;                  % (1/day)
end

% ---------------- K0 CALCULATION -----------

  K0 = exp( ...                                                             % Henrys law constant (mol/kg/atm) or (µmol/kg*ppmv)
       (9345.17/T_K)-60.2409+23.3585*log(T_K/100.) ...
       +S*( 0.023517-0.00023656*T_K  +  ...
       0.00000047036*T_K*T_K )); 

% ----------------ACTIVITY COEFFICIENT CALCULATION-------------%  

% flag options: 
% 'EDH' -> Extened Debye-Hückel Equation
% 'Davies' -> Davies Equation
% 'Ideal' -> Ideal Solution
ac_flag = 'EDH';                                                            % flag to determine how activity coefficients are calculated

lambda_Ca = lambda(S,6,2,ac_flag);                                          % Ca activity coefficient (dimensionless)
lambda_CO3 = lambda(S,4.5,2,ac_flag);                                       % CO3 activity coefficient (dimensionless)
  
% ---------------- PARAMAETER CALCULATION ------------

% ----------------------- FLUX IN --------------------

F_in_DIC    = DIC_in.*Flow_rate.*1E3.*1E-6;                                 % Influent DIC (mol/day) 
F_in_ALK    = ALK_in.*Flow_rate.*1E3.*1E-6;                                 % Influent Alkalinity (mol/day)
F_in_CaCO3  = CaCO3_in*(10^-3)*(1/100.0869)*(10^3)*Flow_rate;               % Influent CaCO3 (mol/day)
F_in_Ca     = Ca_in.*Flow_rate.*1E3.*1E-6;                                  % Influent Dissolved Ca (mol/day)

% ------- Other parameters ------ 

V1   = Time_residence.*Flow_rate;                                           % Volume (m3)
A1   = V1./Depth_average;                                                   % Area (m2)

% -----------------------------------
% ------- SIMULATION  ---------------
% -----------------------------------

t0=0;                                                                       % Initial time (day)

% -------------- INITIAL CONDITIONS -----------------------

% -------------- CO2g ---------------

DIC10    = DIC_in.*1E3.*1E-6.*V1;                                           % Initial DIC (mol)

% -------------- TA ----------------

TA10     = (ALK_in.*1E3.*1E-6.*V1);                                         % Initial Alkalinity (mol)

% ----------- CaCO3 ----------------

CaCO310  = 0;                                                               % Initial CaCO3 (mol)

% -------------- Ca ----------------

Ca0      = (Ca_in.*1E3.*1E-6.*V1);                                          % Initial Ca (mol)

% ----------------------------------

W0 = [DIC10 TA10 CaCO310 Ca0]';                                             % Vector of initial conditions

DIC1conc    = @(W)  W(1)/V1;                                                % DIC concentration (mol / m3)
TA1conc     = @(W)  W(2)/V1;                                                % TA  concentration (mol / m3)
CaCO31conc  = @(W)  W(3);                                                   % CaCO3  quantity (mol)
Ca1conc     = @(W)  W(4)/V1;                                                % Ca  concentration (mol / m3)

% -----------------------------------------------------
% ----------------------- FLUX OUT --------------------
% -----------------------------------------------------

% ----------- DIC --------------

F_out_DIC      = @(W,t)   DIC1conc(W).*Flow_rate;                           % Effluent DIC (mol/day)

% ----------- TA ---------------

F_out_ALK      = @(W,t)   TA1conc(W).*Flow_rate;                            % Effluent Alkalinity (mol/day)

% ----------- Ca ---------------

F_out_Ca       = @(W,t)   Ca1conc(W).*Flow_rate;                            % Effluent Ca (mol/day)

% ----------- CaCO3 ---------------

F_out_CaCO3    = @(W,t)   CaCO31conc(W).*(1/V1).*(ratio_WAS_solids)*(Flow_rate)*(frac_WAS_flow); % Waste sludge CaCO3 (mol/day)

% -----------------------------------------------------
% -------------------- GAS EXCHANGE -------------------
% -----------------------------------------------------

H2CO3_11 = @(W) h2co3_1(TA1conc(W)*(10^6)*(10^-3),DIC1conc(W)*(10^6)*(10^-3),T,S); % H2CO3 Concentration (μmol/L)

% ------- CO2 gas exchange parameters -------

Modern_CO2  = 420;                                                          % Atmospheric CO2 Concentration (ppmv)
kwCO2       = 0.5*24;                                                       % Piston velocity for CO2 (m/day)

F_gas_exchange  = @(W,t) alpha_ebull.*kwCO2.*((H2CO3_11(W)*(10^-6)*(10^3))-(Modern_CO2.*K0.*1E3.*1E-6)).*A1;   % CO2 gas exchange (mol/day)


% ----------------------------------------------------
% ------ Carbonate dissolution/precipitation --------- 
% ----------------------------------------------------

Calcium          = @(W,t)  Ca1conc(W).*1E-3.*1E6;                           % Calcium concentration (mol / L)
Ksp_ca           = ksp_calcite_1(T,S);                                      % Solubility product for calcite (mol2 / m2)

CO3_1            = @(W,t)  co3_1(TA1conc(W).*1E-3.*1E6,DIC1conc(W).*1E-3.*1E6,T,S); % Carbonate ion concentration (μmol / L)
pH_11            = @(W,t)  pH_1(TA1conc(W).*1E-3.*1E6,DIC1conc(W).*1E-3.*1E6,T,S);  % pH (dimensionless)

TDF_diss         = @(W,t)  1-((lambda_Ca*lambda_CO3.*Calcium(W,t).* CO3_1(W,t))./(Ksp_ca)); % Thermodyanmic driving force for dissolution (1-omega) (dimensionless)
TDF_prec         = @(W,t)  ((lambda_Ca*lambda_CO3.*Calcium(W,t).* CO3_1(W,t))./(Ksp_ca))-1; % Thermodyanmic driving force for precipitation (omega-1) (dimensionless)

k_calcite        = 0.0001;                                                  % rate constant for calcite precipitation (1 / day)

heaviside_diss   = @(W,t)  heaviside(TDF_diss(W,t));                        
heaviside_prec   = @(W,t)  heaviside(TDF_prec(W,t));

F_carbonate_prec    = @(W,t)  heaviside_prec(W,t).*-k_calcite.*(TDF_prec(W,t)^n_power_prec).*V1;   % Calcite precipitation rate (mol / day)
F_carbonate_diss    = @(W,t)  heaviside_diss(W,t).*k_calcite_dis.*CaCO31conc(W).*(TDF_diss(W,t)^n_power_diss).*10.^(-n_all.*pH_11(W,t)); % Calcite dissolution rate (mol / day)          

F_carbonate         = @(W,t)  F_carbonate_prec(W,t) +  F_carbonate_diss(W,t); % overall rate of calcite dissolution + precipitation (mol / day)

% ---------------------------------------------------------
% ------------- Denitrification & Nitrification -----------
% ---------------------------------------------------------

NO3_conc = 2.*(1000/62);                                                    % Nitrate Concentration (μM)
NH3_conc = 4.*(1000/17);                                                    % Ammonium Concentration (μM)                                  
O2_conc = 220;                                                              % Oxygen Concentration (μM)

F_nitrif = k_nitrif.*NH3_conc.*O2_conc.*V1.*1E3.*1E-6.*(1/365);             % Nitirification alkalinity removal rate (mol/day)
F_Denitr = R_carbon.*(NO3_conc./(NO3_conc + Dni_half_sat));                 % Denitirification alkalinity addition rate (mol/day)

F_nitrif_a = @(W,t) F_nitrif*(1/(1+exp(-3.5*(pH_11(W,t)-6))));              % Apply logistic pH adjustment to nitrification rate

% ---------------------------------------------------------
% ----------------- SET UP THE ODEs -----------------------
% ---------------------------------------------------------

% ----------- DIC ------------------

dDIC1dt    = @(t,W)  F_in_DIC  - F_out_DIC(W,t) + R_carbon - F_gas_exchange(W,t) + F_carbonate(W,t);  % Rate of change of DIC concentration WRT time (mol / day)

% ----------- TA --------------------

dTA1dt     = @(t,W)  F_in_ALK  - F_out_ALK(W,t) + 2.*F_carbonate(W,t) + F_Denitr - F_nitrif_a(W,t);   % Rate of change of ALK concentration WRT time (mol / day)

% ----------- CaCO3 --------------------

dCaCO31dt  = @(t,W)  F_in_CaCO3 - F_out_CaCO3(W,t) - F_carbonate(W,t) ;                               % Rate of change of calcite concentration WRT time (mol / day)

% -------------- Ca --------------------

dCa1dt      = @(t,W)  F_in_Ca - F_out_Ca(W,t) + F_carbonate(W,t);                                     % Rate of change of calcium concentration WRT time (mol / day) 

% ------------------------------------------------------------------
% -------------- Call the numerical solver  ------------------------
% ------------------------------------------------------------------

% define the vector-valued function of ODEs

dWdt    = @(t,W) [dDIC1dt(t,W) dTA1dt(t,W) dCaCO31dt(t,W) dCa1dt(t,W)]' ;   

% if statement to determine how the model is run: steady state or time
% integraded

if strcmp(run_flag,'steady_state')                                         

    % define the steady state function, here for t = 0

    dWdt_ss = @(W) dWdt(0,W);                       
    
    % set options for fsolve

    opts = optimoptions('fsolve','MaxFunctionEvaluations',10^6,'MaxIterations',2000,"Display","none");     % Extra fsolve options: 'ScaleProblem','jacobian','StepTolerance',10^-10,'FunctionTolerance',10^-10,'MaxFunctionEvaluations',10^6,'MaxIterations',5000,'OptimalityTolerance',10^-20,'Algorithm','levenberg-marquardt','FiniteDifferenceType','central','FiniteDifferenceStepSize',10^-10

    % run fsolve to find steady state

    [W_ss,fval,exitflag,output_message] = fsolve(dWdt_ss,W0,opts);             


elseif strcmp(run_flag,'time_integrate')        

    % set options for the numerical integration

    options = odeset('MaxStep',maxstep); 
    
    % call the stiff ODE solver

    [TT,Y]   = ode15s(dWdt,[t0 tmax],W0,options);                              

end                                         

% ----------------------------------------------------------
% --------------- END OF THE CALCULATION  ------------------
% ----------------------------------------------------------

% if statement to determine if processing time-integrated or steady-state
% results

if  strcmp(run_flag,'time_integrate') 

    DICc   = (Y(:,1).*1E6.*1E-3)./V1;                                           % Effluent DIC (umol / L)
    TA1c   = (Y(:,2).*1E6.*1E-3)./V1;                                           % Effluent Alkalinity (ueq / L)
    CaCO3c = Y(:,3);                                                            % Effluent CaCO3 (mol)      
    Ca1c   = (Y(:,4).*1E6.*1E-3)./V1;                                           % Effluent Dissolved Calcium (mol / L)
    
    % pre-allocate arrays for post-calculations on model results
    
    wCO3 = zeros(size(TT,1),1);
    pH = zeros(size(TT,1),1);
    OMEGA = zeros(size(TT,1),1);
    wF_carbonate = zeros(size(TT,1),1);
    H2CO3_SOLUTION = zeros(size(TT,1),1);
    
    % post calculations using model results
    
        for z=1:size(TT,1)
           wCO3(z)         = CO3_1(Y(z,:),TT(z));
           pH(z)           = pH_1(TA1c(z),DICc(z),T,S);
           OMEGA(z)        = ((Y(z,4)*(10^6)*(1/V1)*(10^-3))*lambda_Ca*lambda_CO3*wCO3(z))./(Ksp_ca);
           wF_carbonate(z) = F_carbonate(Y(z,:),TT(z));
           H2CO3_SOLUTION(z) = H2CO3_11(Y(z,:));
        end
    
    % ---------------------------------------------------------
    % -------------- PLOT RESULTS  ----------------------------
    % ---------------------------------------------------------
    
    linewidth = 2;
    Fontsize = 12;
    
    tiledlayout(4,2)
    
    % DIC
    
    nexttile;
    plot(TT,DICc,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('[DIC] (\muM)' );
    % xlabel ('time (day)');
    set(gca,'XScale','Log')
    % ylim([0.5.*min(DICc) 1.5.*max(DICc)]);
    ylim([0 6000])
    ylim([min(DICc)*0.5 max(DICc)*1.2])
    
    % ALK
    
    nexttile;
    plot(TT,TA1c,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('[ALK] (\muM)' );
    % xlabel ('time (day)');
    set(gca,'XScale','Log')
    % ylim([0.5.*min(TA1c) 1.5.*max(TA1c)]);
    ylim([min(TA1c)*0.5 max(TA1c)*1.2])
    
    % pH
    
    nexttile;
    plot(TT,pH','lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('pH' );
    set(gca,'XScale','Log')
    % xlabel ('time (day)');
    ylim([6 8])
    
    
    % Omega 
    
    nexttile;
    plot(TT,OMEGA,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('\Omega Calcite' );
    set(gca,'XScale','Log')
    % xlabel ('time (day)');
    
    % Solid fraction 
    
    nexttile;
    plot(TT,100.*(CaCO3c./V1),'lineWidth',linewidth);
    % plot(T,wF_carbonate,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('[CaCO_3] (mg/l)' );
    set(gca,'XScale','Log')
    % xlabel ('time (day)');
    
    % Calcium
    
    nexttile;
    plot(TT,Ca1c,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('[Ca] (\muM)' );
    xlabel ('time (days)');
    % ylim([0.5.*min(Ca1c) 1.5.*max(Ca1c)]);
    set(gca,'XScale','Log')
    ylim([min(Ca1c)*0.5 max(Ca1c)*1.2])
    %
    
    %CO2
    
    nexttile;
    plot(TT,H2CO3_SOLUTION./K0,'lineWidth',linewidth);
    set(gca, 'FontName', 'Arial','FontSize', Fontsize,'lineWidth',linewidth);
    ylabel ('[CO2] (\muatm)' );
    xlabel ('time (days)');
    set(gca,'XScale','Log')
    ylim([min(H2CO3_SOLUTION./K0)*0.5 max(H2CO3_SOLUTION./K0)*1.2])
    
    x0=10;
    y0=10;
    width=1200;
    height=800;
    set(gcf,'position',[x0,y0,width,height])

elseif strcmp(run_flag,'steady_state') 

    CDR_Rate = ((W_ss(4)*(1/V1)) - (Ca_in*(10^-6)*(10^3)))*Flow_rate*44.01*(10^-6);
    DIC_ss = W_ss(1)*(1/V1)*(1/1000);
    ALK_ss = W_ss(2)*(1/V1)*(1/1000);
    pH_effluent = pH_1(ALK_ss*(10^6),DIC_ss*(10^6),T,S);

    disp(['Steady-State CDR Rate: ',num2str(CDR_Rate),' (Tonne day-1)'])
    disp(['Steady-State Effluent pH: ',num2str(pH_effluent)])

end