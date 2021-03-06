clear
clc

% ========================== Control Panel ===========================

% ---------- Time Settings ----------

% Time slice
ts = 1/1000; % 10^-3 s

% Total operating time in s
secs = 100;

% ---------- Slag Settings -----------

% Temperature of slag in K
T_slg = 300;

% Slag mass addition rate in kg/s
slg_add = 3;

% Slag mass fraction
MX_CaO_slg = 0.573;
MX_MgO_slg = 0.415;
MX_SiO2_slg = 0.007;
MX_Al2O3_slg = 0.005;

% --------- Reactor Geometry ---------
r_eafout = 3.3;
r_eafin = 2.45;
r_hole = 1.7;
r_electrode = 0.3;
h_eafup = 2.9;
h_eaflow = 1.0;
h_electrode = 1.0;
d1 = 0.30;
d2 = 0.45;

% ---------- Other Settings ----------

% Carbon Injection Rate (kg/s)
C_inj = 0;

% Manganese Injection Rate (kg/s)
Mn_inj = 0;

% Oxygen Lance Rate (kg/s)
O2_lance = 3.2;

% O2 for post combustion (kg/s)
O2_post = 0;

% Power of arc (kW)
P_arc = 48900;

% EAF mass capacity (kg)
m_EAF = 105000;

% Cooling water flowrate (mol/s)
phi1 = 80/0.018;
phi2 = 150/0.018;

% ======================= Initial Parameters =========================

% ------------- Initial mass (kg) --------------

% Solid metal initial mass
m_Fe_sSc = 7764;
m_C_sSc = 32;
m_Cr_sSc = 16;
m_Mn_sSc = 48;
m_P_sSc = 4;
m_Si_sSc = 48;
m_comb_sSc = 88;

m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + ...
    m_Si_sSc + m_comb_sSc;

% Liquid metal initial mass
m_Fe_lSc = 31442;
m_C_lSc = 112;
m_Cr_lSc = 51.2;
m_Mn_lSc  = 186;
m_P_lSc = 22.4;
m_Si_lSc = 186;

m_lSc = m_Fe_lSc + m_C_lSc + m_Cr_lSc + m_Mn_lSc + m_P_lSc + m_Si_lSc;

% Solid slag initial mass
m_CaO_sSl = 113;
m_MgO_sSl = 109;
m_SiO2_sSl = 2;
m_Al2O3_sSl = 3.6;

% Liquid slag initial mass
m_SiO2_lSl = 307;
m_Al2O3_lSl = 153;
m_CaO_lSl = 313;
m_MgO_lSl = 395;
m_MnO_lSl = 51;
m_P2O5_lSl = 11;
m_Cr2O3_lSl = 2.19;
m_FeO_lSl = 310;

m_lSl = m_SiO2_lSl + m_Al2O3_lSl + m_CaO_lSl + m_MgO_lSl + m_MnO_lSl + ...
    m_P2O5_lSl + m_Cr2O3_lSl + m_FeO_lSl;

% Gas initial mass
m_H2O = 18;
m_O2 = 41;
m_CO = 80;
m_CO2 = 19;

% Initial mass of injected carbon
m_CL = 0.0671;

% ------------- Initial temp. (K) --------------

T_sSc = 1100;
T_lSc = 1800;
T_sSl = 1100;
T_lSl = 1800;
T_gas = 2000;
T_wall = 328;
T_roof = 344;

% -------------- Initial Geometry --------------

% Densities (kg/m3)
rho_sSc = 900; % Kurz and Fisher 2005
rho_lSc = 7000;
rho_lSl = 3500; % Self-compacting concrete: materials properties and applications Siddique (2020) Table 10.1
rho = 7000;

% Cross-sectional areas
A_bath = pi * r_eafin^2;
A_eaf = pi * r_eafout^2;

% height of liquid metal
h_lSc = (m_lSc / rho_lSc) / A_bath;

% height of liquid slag
h_lSl = (m_lSl / rho_lSl) / A_bath;

% height of solid metal
% V_free = (A_bath*h_eaflow) - (m_lSc/rho_lSc) - (m_lSl/rho_lSl);
% if V_free >= m_sSc/rho_sSc
%     h_sSc1 = 0;
%     h_sSc2 = (m_sSc/rho_sSc) / A_bath;
% else
%     h_sSc1 = ((m_sSc/rho_sSc) - V_free) / A_eaf;
%     h_sSc2 = V_free / A_bath;
% end
% 

h_sSc2 = 0;
d_conein = r_eafin*2;
d_coneout = r_eafout*2;

V_sSc = m_sSc / rho_sSc;
h_cone = V_sSc / (pi*r_eafout^2-(1/3)*pi*(r_eafin^2 + r_eafin*r_eafout + r_eafout^2));
h_sSc1 = h_cone;

% Arc height
h_arc = h_eafup - h_electrode - (h_sSc2 - h_cone);

% height of wall
h_wall = h_eafup + h_eaflow - h_sSc1 - h_lSc - h_lSl;

% Areas of roof and wall
A1 = (pi * r_eafout^2) - (pi * r_hole^2); % roof
A2 = 2 * pi * r_eafout * h_wall; % wall
A4 = pi * r_eafin^2;

% ----------------- Others ---------------------

% Initial Pressure
p_gas = 1.2; % atm
rp = 0; % Relative pressure

% Exposure constant
K_sSclSc = 0.5*tanh(5*(h_lSc-h_sSc1-h_sSc2+h_cone)) + 0.5;

% Initial Radiosity
J_roof = 0;
J_wall = 0;
J_sSc = 0;
J_lSc = 0;

% =========================== Constants ==============================

% ------------ Molar Mass (kg/mol) -------------
M_Fe = 0.055845;
M_C = 0.0120107;
M_Si = 0.0280855; 
M_FeO = 0.071844;
M_SiO2 = 0.06008; 
M_O2 = 0.032;
M_CO = 0.02801;
M_CO2 = 0.04401;
M_CH4 = 0.01604;
M_Cr = 0.051996;
M_Mn = 0.054938;
M_P = 0.030976;
M_MnO = 0.0709374;
M_Cr2O3 = 0.15199;
M_P2O5 = 0.283886 / 2;
M_CaO = 0.05608;
M_MgO = 0.04304;
M_C9H20 = 0.1282;
M_gas = 0.035;
M_Al2O3 = 0.10196;
M_Al = 0.02698;
M_H2O = 0.01802;
M_sSl = 0.0509;
M_lSl = 0.0712;

% ------------ Reaction rate constant -------------

kd_CL = 15; % s-1
kd_CD = 35; % kg/s
kd_C1 = 60;
kd_C2 = 55;
kd_Mn1 = 20;
kd_Mn2 = 10;
kd_Mn = 75;
kd_Si1 = 144;
kd_Si2 = 250;
kd_Cr1 = 3;
kd_Cr2 = 1;
kd_P = 35;
kd_comb = 0.1; %s-1

% ----------- Thermodynamic Properties ------------

% Heat capacities (kJ/mol K)
Cp_C = 0.02092;  % at 1800K (from NIST)
Cp_H2O = 0.075;
Cp_sSc = 0.039;
Cp_lSc = 0.047;
Cp_sSl = 0.025;
Cp_lSl = 0.047;
Cp_gas = 0.030;
Cp_roof = 0.65; % kJ/kg K
Cp_wall = 0.96; % kJ/kg K

% Latent heat of fusion (kJ/mol)
lambda_C = 117;
lambda_sSc = 15.4;
lambda_sSl = 12.66;

% Enthalpies of formation (kJ/mol)
dH_FeO = -243;
dH_FeS = 0;
dH_CO = -117;
dH_CO2 = -396;
dH_CS = -27;
dH_MnS = -20;
dH_MnO = -385;
dH_SiS = -132;
dH_SiO2 = -946;
dH_SiO2S = -45;
dH_CrS = -42;
dH_Cr2O3 = -1142;
dH_PS = -29;
dH_P2O5 = -2940;
dH_H2O = -247;
dH_CH4 = -91;
dH_C9H20 = -228.3;

% ------------ Conduction Constants ---------------

% Heat transfer coeffs kW/m^2K
K_therm1 = 0.2;
K_therm2 = 0.2;
K_therm3 = 0.05;
K_therm4 = 57.5; % kW/K
K_therm5 = 0.08;
K_therm6 = 0.08;
K_therm7 = 22.5; % kW/K
K_therm8 = 22.5; % kW/K

% Area coeffs m^2/kg
K_area1 = 0.008;
K_area2 = 0.12;
K_area3 = 0.12;
K_area5 = 0.12;
K_area6 = 0.12;

% Thermal conductance coeff. kW/K
K_water1 = 12;
K_water2 = 20;
K_water3 = 10;
K_water4 = 5;
K_water5 = 0.05;

% ------------ Radiation Constants ---------------

% Stefan Boltzmann Constant
sig = 5.67e-08; % W.m-2.K-4s-1

% Emissivity
ep1 = 0.85; % Roof
ep2 = 0.85; % Wall
ep3 = 0.80; % sSc
ep4 = 0.40; % lSc

% ------------ Constant Temperatures --------------

T_air = 298;
T_melt = 1809;
T_water = 298;

% ------------------- Others ----------------------

% Distribution of lanced oxygen
K_O2CO = 0.05;
K_O2CO2 = 0.15;
K_O2Cr2O3 = 0.015;
K_O2FeO = 0.75;
K_O2SiO2 = 0.035;

% Vent
hd = 0.65;
k_U = 6.44;
u1 = 17;
u2 = 0.3;

V_gas = 45;

% Electrode
R_tip = 0.02;
R_side = 10;
A_side = 35;
I_arc = 44;

% ------------------------- Arrays for graph ------------------------
gas_temp = zeros(1,secs*5);
sSc_temp = zeros(1,secs*5);
sSl_temp = zeros(1,secs*5);
lSc_temp = zeros(1,secs*5);
lSl_temp = zeros(1,secs*5);
steel_Fe = zeros(1,secs*5);
m_solid = zeros(1,secs*5);
m_liquid = zeros(1,secs*5);
m_solid_slag = zeros(1,secs*5);
m_liquid_slag = zeros(1,secs*5);
rel_pres = zeros(1,secs*5);

for step = 1:secs/ts
    
    % ====================== Mole, Mass Fraction =====================
    
    % ------------ Solid Metal ------------
    
    % Total mass of solid metal
    m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_Si_sSc + ...
        m_comb_sSc;
    
    % Total mole of solid metal
    XM_sSc = (m_Fe_sSc/M_Fe) + (m_C_sSc/M_C) + (m_Cr_sSc/M_Cr) + (m_Mn_sSc/M_Mn) ... 
        + (m_P_sSc/M_P) + (m_Si_sSc/M_Si) + (m_comb_sSc/M_C9H20);
    
    % Mole fractions of compounds in solid metal
    X_Fe_sSc = (m_Fe_sSc/M_Fe) / XM_sSc;
    X_C_sSc = (m_C_sSc/M_C) / XM_sSc;
    X_Cr_sSc = (m_Cr_sSc/M_Cr) / XM_sSc;
    X_Mn_sSc = (m_Mn_sSc/M_Mn) / XM_sSc;
    X_P_sSc = (m_P_sSc/M_P) / XM_sSc;
    X_Si_sSc = (m_Si_sSc/M_Si) / XM_sSc;
    X_comb_sSc = (m_comb_sSc/M_C9H20) / XM_sSc;
    
    % Mass fractions of compounds in solid metal
    MX_Fe_sSc = m_Fe_sSc / m_sSc;
    MX_C_sSc = m_C_sSc / m_sSc;
    MX_Cr_sSc = m_Cr_sSc / m_sSc;
    MX_Mn_sSc = m_Mn_sSc / m_sSc;
    MX_P_sSc = m_P_sSc / m_sSc;
    MX_Si_sSc = m_Si_sSc / m_sSc;
    MX_comb_sSc = m_comb_sSc / m_sSc;
    
    % ----------- Liquid Metal ------------
    
    % Total mass of liquid metal
    m_lSc = m_Fe_lSc + m_C_lSc + m_Cr_lSc + m_Mn_lSc + m_P_lSc + m_Si_lSc;
    
    % Total mole of liquid metal
    XM_lSc = (m_Fe_lSc/M_Fe) + (m_C_lSc/M_C) + (m_Cr_lSc/M_Cr) + (m_Mn_lSc/M_Mn) ...
        + (m_P_lSc/M_P) + (m_Si_lSc/M_Si);
    
    % Mole fractions of compounds in liquid metal
    X_Fe_lSc = (m_Fe_lSc/M_Fe) / XM_lSc;
    X_C_lSc = (m_C_lSc/M_C) / XM_lSc;
    X_Cr_lSc = (m_Cr_lSc/M_Cr) / XM_lSc;
    X_Mn_lSc = (m_Mn_lSc/M_Mn) / XM_lSc;
    X_P_lSc = (m_P_lSc/M_P) / XM_lSc;
    X_Si_lSc = (m_Si_lSc/M_Si) / XM_lSc;
    
    % Mass fractions of compounds in liquid metal
    MX_Fe_lSc = m_Fe_lSc / m_lSc;
    MX_C_lSc = m_C_lSc / m_lSc;
    MX_Cr_lSc = m_Cr_lSc / m_lSc;
    MX_Mn_lSc = m_Mn_lSc / m_lSc;
    MX_P_lSc = m_P_lSc / m_lSc;
    MX_Si_lSc = m_Si_lSc / m_lSc;
    
    % ------------ Solid Slag -------------
    
    % Total mass of solid slag
    m_sSl = m_CaO_sSl + m_MgO_sSl + m_SiO2_sSl + m_Al2O3_sSl;
    
    % Total mole of solid slag
    XM_sSl = (m_CaO_sSl/M_CaO) + (m_MgO_sSl/M_MgO) + (m_SiO2_sSl/M_SiO2) ...
        + (m_Al2O3_sSl/M_Al2O3);
    
    % Mole fractions of compounds in solid slag
    X_CaO_sSl = (m_CaO_sSl/M_CaO) / XM_sSl; 
    X_MgO_sSl = (m_MgO_sSl/M_MgO) / XM_sSl;
    X_SiO2_sSl = (m_SiO2_sSl/M_SiO2) / XM_sSl;
    X_Al2O3_sSl = (m_Al2O3_sSl/M_Al2O3) / XM_sSl;
    
    % Mass fractions of compounds in solid slag
    MX_CaO_sSl = m_CaO_sSl / m_sSl;
    MX_MgO_sSl = m_MgO_sSl / m_sSl;
    MX_SiO2_sSl = m_SiO2_sSl / m_sSl;
    MX_Al2O3_sSl = m_Al2O3_sSl / m_sSl;
    
    % ----------- Liquid Slag -------------
    
    % Total mass of liquid slag
    m_lSl = m_SiO2_lSl + m_Al2O3_lSl + m_CaO_lSl + m_MgO_lSl + m_MnO_lSl ...
        + m_P2O5_lSl + m_Cr2O3_lSl + m_FeO_lSl;
    
    % Total mole of liquid slag
    XM_lSl = (m_SiO2_lSl/M_SiO2) + (m_Al2O3_lSl/M_Al2O3) + (m_CaO_lSl/M_CaO) ...
        + (m_MgO_lSl/M_MgO) + (m_MnO_lSl/M_MnO) + (m_P2O5_lSl/M_P2O5) + (m_Cr2O3_lSl/M_Cr2O3) ...
        + (m_FeO_lSl/M_FeO);
    
    % Mole fractions of compounds in liquid slag
    X_SiO2_lSl = (m_SiO2_lSl/M_SiO2) / XM_lSl;
    X_Al2O3_lSl = (m_Al2O3_lSl/M_Al2O3) / XM_lSl;
    X_CaO_lSl = (m_CaO_lSl/M_CaO) / XM_lSl;
    X_MgO_lSl = (m_MgO_lSl/M_MgO) / XM_lSl;
    X_MnO_lSl = (m_MnO_lSl/M_MnO) / XM_lSl;
    X_P2O5_lSl = (m_P2O5_lSl/M_P2O5) / XM_lSl;
    X_Cr2O3_lSl = (m_Cr2O3_lSl/M_Cr2O3) / XM_lSl;
    X_FeO_lSl = (m_FeO_lSl/M_FeO) / XM_lSl;
    
    % Mass fractions of compounds in solid slag
    MX_SiO2_lSl = m_SiO2_lSl / m_lSl;
    MX_Al2O3_lSl = m_Al2O3_lSl / m_lSl;
    MX_CaO_lSl = m_CaO_lSl / m_lSl;
    MX_MgO_lSl = m_MgO_lSl / m_lSl;
    MX_MnO_lSl = m_MnO_lSl / m_lSl;
    MX_P2O5_lSl = m_P2O5_lSl / m_lSl;
    MX_Cr2O3_lSl = m_Cr2O3_lSl / m_lSl;
    MX_FeO_lSl = m_FeO_lSl / m_lSl;
    
    % --------------- Gas -----------------
    
    % Total mass of gas
    m_gas = 80;
    
    % Total mole of gas
    XM_gas = (m_H2O/M_H2O) + (m_CO/M_CO) + (m_CO2/M_CO2) + (m_O2/M_O2);
    
    % Mole fraction of compounds in gas
    X_H2O = (m_H2O/M_H2O) / XM_gas;
    X_CO = (m_CO/M_CO) / XM_gas;
    X_CO2 = (m_CO2/M_CO2) / XM_gas;
    X_O2 = (m_O2/M_O2) / XM_gas;
    
    % Mass fraction of compounds in gas
    MX_H2O = m_H2O / m_gas;
    MX_CO = m_CO / m_gas;
    MX_CO2 = m_CO2 / m_gas;
    MX_O2 = m_O2 / m_gas;
    
    % ====================== Chemical Reactions ======================
    
    % ------- Injected Carbon Reaction -------
    
    % FeO + C -> Fe + CO
    % Decarburization reaction in mol/s
    if isnan(MX_FeO_lSl)
        r_FeO_CL = 0;
    else
        r_FeO_CL = (m_CL * kd_CL * MX_FeO_lSl) / M_C;
    end
    
    % ------ Dissolved C decarburization -------
    
    % FeO + C -> Fe + CO
    
    Xeq_C = 4.9e-04 / X_FeO_lSl;
    if isnan(Xeq_C)
        Xeq_C = 0;
    end
    % Rate of decarburization of dissolved C
    r_FeO_CD = (kd_CD * (X_C_lSc - Xeq_C)) / M_C; % Logar 2012
    
    % ---------- Carbon Oxidation -----------
    
    % To carbon monoxide
    % C + 1/2 O2 -> CO
    r_C_hO2 = (kd_C1 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO) / M_C; % Logar 2012
    
    % To carbon dioxide
    % C + O2 -> CO2
    r_C_O2 = (kd_C2 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO2) / M_C; % Logar 2012
    
    % -------- MnO decarburization ----------
    
    % MnO + C -> Mn + CO
    
    p_CO = p_gas * X_CO;
    
    % Equilibrium constant
    kX_Mn1 = 6.4 * p_CO * X_MnO_lSl;
    
    % Equilibrium mole fraction
    Xeq_MnO1 = X_Mn_lSc / kX_Mn1; % Logar 2012
    
    if isnan(Xeq_MnO1)
        Xeq_MnO1 = 0;
    end
    
    % Rate of reaction
    r_MnO_C = (kd_Mn1 * (X_MnO_lSl - Xeq_MnO1)) / M_MnO;
    
    % ----------- Desiliconization ----------
    
    % 2FeO + Si -> 2Fe + SiO2
    
    % Basicity
    B3 = (MX_CaO_lSl * MX_MgO_lSl) / (MX_SiO2_lSl * MX_Al2O3_lSl);
    
    % Activity of SiO2 based on basicity
    a_SiO2_bas = exp((6728/T_lSl) - (0.920*B3 + 6.994)); % Meraikib 1995
    
    % Equilibrium constant
    K_Si = 10 ^ (30410/T_lSc - 11.59); % Turkdogan 1996
    
    % Oxygen Solubility
    O_sol = (10 ^ (-6380/T_lSc + 2.765)); % Turkdogan 1996
    
    % Equilibrium fraction
    MXeq_Si = a_SiO2_bas / (K_Si * O_sol^2); % Turkogan 1996
    Xeq_Si = MXeq_Si * (M_Fe/(M_Si*100));
    % Xeq_Si = ((MXeq_Si * m_lSc) / M_Si) / XM_lSc; 
    % Xeq_Si = 8.08e-08* ((m_lSl*M_FeO)/(m_FeO_lSl*M_lSl) + (m_SiO2_lSl*M_FeO) ...
    %     /(m_FeO_lSl*M_SiO2) + 1)^2; % Bekkar 1999
    
    if isnan(Xeq_Si)
        Xeq_Si = 0;
    end

    % Rate of reaction
    r_2FeO_Si = (kd_Si1 * (X_Si_lSc - Xeq_Si)) / M_Si;
    
    % --------- Silicon Oxidation ---------
    
    % Si + O2 -> SiO2
    
    r_Si_O2 = (kd_Si2 * (X_Si_lSc - Xeq_Si) * O2_lance * K_O2SiO2) / M_Si;
    
    % ------- Si reaction with MnO --------
    
    % 2MnO + Si -> 2Mn + SiO2
    
    % Reaction equilibrium constant
    kX_Mn2 = 10^(2.8*((X_CaO_lSl+X_MgO_lSl)/X_SiO2_lSl)/-1.16) * ((M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2)); % Logar 2012
    
    % Equilibrium MnO mole fraction
    Xeq_MnO2 = 0.01; % Logar 2012
    
    if isnan(Xeq_MnO2)
        Xeq_MnO2 = 0;
    end
    
    % Rate of reaction
    r_2MnO_Si = (kd_Mn2 * (X_MnO_lSl - Xeq_MnO2)) / M_MnO;
    
    % -------- Mn reaction with FeO --------
    
    % Mn + FeO -> MnO + Fe
    
    % Equilibrium constant
    K_FeMn = 1.8; % Turkdogan 1996 (Given that B = 2.5-4.0 and 1600 - 1650 C)
    kX_Mn = K_FeMn * (M_FeO * M_Mn * 100) / (M_MnO * M_Fe);
    
    % Equilibrium fraction
    Xeq_Mn = X_MnO_lSl / (X_FeO_lSl * kX_Mn);
    
    if isnan(Xeq_Mn)
        Xeq_Mn = 0;
    end
    
    % Rate of reaction
    r_FeO_Mn = (kd_Mn * (X_Mn_lSc - Xeq_Mn)) / M_Mn;
    
    % ------- Cr reaction with FeO ---------
    
    % 3FeO + 2Cr -> 3Fe + Cr2O3
    
    % Equilibrium constants
    K_FeCr = 0.3; % Trukdogan
    kX_Cr = K_FeCr * (M_Cr*M_FeO*100) / (M_Cr2O3*M_Fe);
    
    % Equilibrium mole fraction
    Xeq_Cr = X_Cr2O3_lSl / (X_FeO_lSl*kX_Cr);
    
    if isnan(Xeq_Cr)
        Xeq_Cr = 0;
    end
    
    % Rate of reaction
    r_3FeO_2Cr = (kd_Cr1 * (X_Cr_lSc - Xeq_Cr)) / M_Cr;
    
    % --------- Chromium Oxidation ---------
    
    % 2Cr + 3/2O2 -> Cr2O3
    
    r_2Cr_3hO2 = (2*kd_Cr2 * (X_Cr_lSc - Xeq_Cr) * O2_lance * K_O2Cr2O3) / M_Cr;
    
    % ----- Phosphorus reaction with FeO -----
    
    % 5FeO + 2P -> 5Fe + P2O5
    
    Xeq_P = 0.001;
    r_5FeO_2P = (kd_P * (X_P_lSc - Xeq_P)) / M_P;
    if r_5FeO_2P < 0
        r_5FeO_2P = 0;
    end
    
    % ------------- Fe Oxidation -------------
    
    % Fe + 1/2O2 -> FeO
    
    r_Fe_hO2 = (2 * O2_lance * K_O2FeO) / M_O2;
    
    % ------------- Combustion ---------------
    
    % C9H20 + 14O2 -> 9CO2 + 10H2O
    
    r_comb = (kd_comb * m_comb_sSc * (T_sSc/T_melt)) / M_C9H20;
    
    % ----------- Post Combustion ------------
    
    % CO + 1/2O2 -> CO2
    
    K_mCO = 0.9;
    
    r_post = (O2_post * K_mCO) / M_O2;
    
    % --------- Electrode Oxidation ----------
    
    % C + O2 -> CO2
    
    dm_el = 3*(R_tip * (I_arc^2/3600) + R_side * (A_side/3600));
    
    % =================== Reaction Heat Transfer ====================
    
    % ------------ Dynamic heat capacity calculation ----------
    % Units in kJ/mol K
    
    t = (T_gas-298)/1000;
    if T_gas < 500
        CpdT_H2O = -203.6060*t + (1523.290/2)*t^2 + (-3196.413/3)*t^3 + (2474.455/4)*t^4 ...
            - (3.855326)/t;
    elseif T_gas < 1700
        t = (500-298)/1000;
        CpdT_H2O1 = -203.6060*t + (1523.290/2)*t^2 + (-3196.413/3)*t^3 + (2474.455/4)*t^4 ...
            - (3.855326)/t;
        t = (T_gas-500)/1000;
        CpdT_H2O2 = 30.0920*t + (6.8532514/2)*t^2 + (6.793435/3)*t^3 + (-2.534480/4)*t^4 ...
            - (0.082139)/t;
        CpdT_H2O = CpdT_H2O1 + CpdT_H2O2;
    else
        t = (500-298)/1000;
        CpdT_H2O1 = -203.6060*t + (1523.290/2)*t^2 + (-3196.413/3)*t^3 + (2474.455/4)*t^4 ...
            - (3.855326)/t;
        t = (1700-500)/1000;
        CpdT_H2O2 = 30.0920*t + (6.8532514/2)*t^2 + (6.793435/3)*t^3 + (-2.534480/4)*t^4 ...
            - (0.082139)/t;
        t = (T_gas-1700)/1000;
        CpdT_H2O3 = 41.96426*t + (8.622053/2)*t^2 + (-1.499780/3)*t^3 + (0.098119/4)*t^4 ...
            - (0.082139)/t;
        CpdT_H2O = CpdT_H2O1 + CpdT_H2O2 + CpdT_H2O3;
    end
    
    t = (T_lSc-298)/1000;
    CpdT_C = 21.17510*t + (-0.812428/2)*t^2 + (0.448537/3)*t^3 + (-0.043256/4)*t^4 ...
            - (-0.013103)/t;
        
    t = (T_gas-298)/1000;
    CpdT_C_gas = 21.17510*t + (-0.812428/2)*t^2 + (0.448537/3)*t^3 + (-0.043256/4)*t^4 ...
            - (-0.013103)/t;
    
    if T_lSc < 1650
        CpdT_FeO = 45.75120*t + (18.78553/2)*t^2 + (-5.952201/3)*t^3 + (0.852779/4)*t^4 ...
            - (-0.081265)/t;
    else
        t = (1650-298)/1000;
        CpdT_FeO1 = 45.75120*t + (18.78553/2)*t^2 + (-5.952201/3)*t^3 + (0.852779/4)*t^4 ...
            - (-0.081265)/t;
        t = (T_lSc-1650)/1000;
        CpdT_FeO2 = 68.1992*t;
        CpdT_FeO = CpdT_FeO1 + CpdT_FeO2;
    end
    
    t = (T_lSc-298)/1000;
    CpdT_Fe = 23.97449*t + (8.367750/2)*t^2 + (0.000277/3)*t^3 + (-0.000086/4)*t^4 ...
            - (-0.000005)/t;
    
    if T_gas < 700
        t = (T_gas-298)/1000;
        CpdT_O2_gas = 31.32234*t + (-20.23532/2)*t^2 + (57.86644/3)*t^3 + (-36.50624/4)*t^4 ...
            - (-0.007374)/t;
    else
        t = (700-298)/1000;
        CpdT_O21 = 31.32234*t + (-20.23532/2)*t^2 + (57.86644/3)*t^3 + (-36.50624/4)*t^4 ...
            - (-0.007374)/t;
        t = (T_gas-700)/1000;
        CpdT_O22 = 30.03235*t + (8.772972/2)*t^2 + (-3.988133/3)*t^3 + (0.788313/4)*t^4 ...
            - (-0.741599)/t;
        CpdT_O2_gas = CpdT_O21 + CpdT_O22;
    end
    
    if T_lSc < 700
        t = (T_lSc-298)/1000;
        CpdT_O2_gas = 31.32234*t + (-20.23532/2)*t^2 + (57.86644/3)*t^3 + (-36.50624/4)*t^4 ...
            - (-0.007374)/t;
    else
        t = (700-298)/1000;
        CpdT_O21 = 31.32234*t + (-20.23532/2)*t^2 + (57.86644/3)*t^3 + (-36.50624/4)*t^4 ...
            - (-0.007374)/t;
        t = (T_lSc-700)/1000;
        CpdT_O22 = 30.03235*t + (8.772972/2)*t^2 + (-3.988133/3)*t^3 + (0.788313/4)*t^4 ...
            - (-0.741599)/t;
        CpdT_O2_lSc = CpdT_O21 + CpdT_O22;
    end
    
    if T_gas < 1300
        t = (T_gas-298)/1000;
        CpdT_CO_gas = 25.56759*t + (6.096130/2)*t^2 + (4.054656/3)*t^3 + (-2.671301/4)*t^4 ...
            - (-2.671301)/t;
    else
        t = (1300-298)/1000;
        CpdT_CO1 = 25.56759*t + (6.096130/2)*t^2 + (4.054656/3)*t^3 + (-2.671301/4)*t^4 ...
            - (-2.671301)/t;
        t = (T_gas-1300)/1000;
        CpdT_CO2 = 35.15070*t + (1.300095/2)*t^2 + (-0.205921/3)*t^3 + (0.013550/4)*t^4 ...
            - (-3.282780)/t;
        CpdT_CO_gas = CpdT_CO1 + CpdT_CO2;
    end
    
    if T_lSc < 1300
        t = (T_lSc-298)/1000;
        CpdT_CO_lSc = 25.56759*t + (6.096130/2)*t^2 + (4.054656/3)*t^3 + (-2.671301/4)*t^4 ...
            - (-2.671301)/t;
    else
        t = (1300-298)/1000;
        CpdT_CO1 = 25.56759*t + (6.096130/2)*t^2 + (4.054656/3)*t^3 + (-2.671301/4)*t^4 ...
            - (-2.671301)/t;
        t = (T_lSc-1300)/1000;
        CpdT_CO2 = 35.15070*t + (1.300095/2)*t^2 + (-0.205921/3)*t^3 + (0.013550/4)*t^4 ...
            - (-3.282780)/t;
        CpdT_CO_lSc = CpdT_CO1 + CpdT_CO2;
    end
    
    CpdT_MnO = 0.04869 * (T_lSc - 298);
    CpdT_Mn = 0.040 * (T_lSc - 298);
    
    if T_lSc < 847
        t = (T_lSc - 298)/1000;
        CpdT_SiO2 = -6.076591*t + (251.6755/2)*t^2 + (-324.7964/3)*t^3 + (168.5604/4)*t^4 ...
            - (0.002548)/t;
    else
        t = (847 - 298)/1000;
        CpdT_SiO21 = -6.076591*t + (251.6755/2)*t^2 + (-324.7964/3)*t^3 + (168.5604/4)*t^4 ...
            - (0.002548)/t;
        t = (T_lSc - 847)/1000;
        CpdT_SiO22 = 58.75340*t + (10.27925/2)*t^2 + (-0.131384/3)*t^3 + (0.025210/4)*t^4 ...
            - (0.025601)/t;
        CpdT_SiO2 = CpdT_SiO21 + CpdT_SiO22;
    end
    
    t = (T_lSc - 298)/1000;
    CpdT_Si = 22.81719*t + (3.899510/2)*t^2 + (-0.082885/3)*t^3 + (0.042111/4)*t^4 ...
            - (-0.354063)/t;
    
    CpdT_Cr2O3 = 124.6550*t + (-0.337045/2)*t^2 + (5.705010/3)*t^3 + (-1.053470/4)*t^4 ...
            - (-2.030501)/t;
    
    if T_lSc < 600
        t = (T_lSc - 298)/1000;
        CpdT_Cr = 7.489737*t + (71.50498/2)*t^2 + (-91.67562/3)*t^3 + (46.04450/4)*t^4 ...
            - (0.138157)/t;
    else
        t = (600 - 298)/1000;
        CpdT_Cr1 = 7.489737*t + (71.50498/2)*t^2 + (-91.67562/3)*t^3 + (46.04450/4)*t^4 ...
            - (0.138157)/t;
        t = (T_lSc - 600)/1000;
        CpdT_Cr2 = 18.46508*t + (5.477986/2)*t^2 + (7.904329/3)*t^3 + (-1.147848/4)*t^4 ...
            - (1.265791)/t;
        CpdT_Cr = CpdT_Cr1 + CpdT_Cr2;
    end
    
    CpdT_P2O5 = 0.143 * (T_lSc - 298);
    
    if T_lSc < 1180
        CpdT_P = 0.02633 * (T_lSc -298);
    else
        CpdT_P1 = 0.02633 * (1180 -298);
        t = (T_lSc - 1180)/1000;
        CpdT_P2 = 20.44403*t + (1.051745/2)*t^2 + (-1.098514/3)*t^3 + (0.377924/4)*t^4 ...
            - (0.010645)/t;
        CpdT_P = CpdT_P1 + CpdT_P2;
    end
    
    if T_gas < 1200
        t = (T_gas - 298)/1000;
        CpdT_CO2_gas = 24.99735*t + (55.18696/2)*t^2 + (-33.69137/3)*t^3 + (7.948387/4)*t^4 ...
            - (-0.136638)/t;
    else
        t = (1200 - 298)/1000;
        CpdT_CO21 = 24.99735*t + (55.18696/2)*t^2 + (-33.69137/3)*t^3 + (7.948387/4)*t^4 ...
            - (-0.136638)/t;
        t = (T_gas - 1200)/1000;
        CpdT_CO22 = 58.16639*t + (2.720074/2)*t^2 + (-0.492289/3)*t^3 + (0.038844/4)*t^4 ...
            - (-6.447293)/t;
        CpdT_CO2_gas = CpdT_CO21 + CpdT_CO22;
    end
    if T_lSc < 1200
        t = (T_lSc - 298)/1000;
        CpdT_CO2_lSc = 24.99735*t + (55.18696/2)*t^2 + (-33.69137/3)*t^3 + (7.948387/4)*t^4 ...
            - (-0.136638)/t;
    else
        t = (1200 - 298)/1000;
        CpdT_CO21 = 24.99735*t + (55.18696/2)*t^2 + (-33.69137/3)*t^3 + (7.948387/4)*t^4 ...
            - (-0.136638)/t;
        t = (T_lSc - 1200)/1000;
        CpdT_CO22 = 58.16639*t + (2.720074/2)*t^2 + (-0.492289/3)*t^3 + (0.038844/4)*t^4 ...
            - (-6.447293)/t;
        CpdT_CO2_lSc = CpdT_CO21 + CpdT_CO22;
    end
    Cp_CH4 = 0.0586;
    CpdT_C9H20 = 0.40334 * (T_gas - 298);
    
    % ----------------- Heat of reaction -----------------
    
    % a) Fe + 1/2O2 -> FeO
    dH_Ta = r_Fe_hO2 * (dH_FeO + CpdT_FeO - CpdT_Fe - 0.5*CpdT_O2_lSc);

    % b) FeO + C -> Fe + CO
    dH_Tb = (r_FeO_CL + r_FeO_CD) * (dH_CO - dH_CS - dH_FeO + CpdT_Fe + CpdT_CO_lSc - CpdT_C - CpdT_FeO);

    % c) FeO + Mn -> Fe + MnO
    dH_Tc = r_FeO_Mn * (dH_MnO - dH_FeO - dH_MnS + CpdT_Fe + CpdT_MnO - CpdT_FeO - CpdT_Mn);

    % d) 2FeO + Si -> 2Fe + SiO2
    dH_Td = r_2FeO_Si * ((dH_SiO2+dH_SiO2S-2*dH_FeO-dH_SiS) + ...
        (2*CpdT_Fe + CpdT_SiO2 - 2*CpdT_FeO - CpdT_Si));

    % e) 3FeO + 2Cr -> 3Fe + Cr2O3
    dH_Te = r_3FeO_2Cr * (dH_Cr2O3 - 3*dH_FeO - 2*dH_CrS + ...
        3*CpdT_Fe + CpdT_Cr2O3 - 3*CpdT_FeO - 2*CpdT_Cr);

    % f) 5FeO + 2P -> 5Fe + P2O5
    dH_Tf = r_5FeO_2P * (dH_P2O5 - 5*dH_FeO - 2*dH_PS + ...
        5*CpdT_Fe + CpdT_P2O5 - 5*CpdT_FeO - 2*CpdT_P);

    % g) C + 1/2O2 -> CO
    dH_Tg = r_C_hO2 * ((dH_CO-dH_CS) + CpdT_CO_lSc - CpdT_C - CpdT_O2_lSc);

    % h) CO + 1/2O2 -> CO2
    dH_Th = r_post * ((dH_CO2-dH_CO) + CpdT_CO2_gas - CpdT_CO_gas - 0.5*CpdT_O2_gas);

    % i) C + O2 -> CO2
    dH_Ti = r_C_O2 * ((dH_CO2-dH_CS) + CpdT_CO2_lSc - CpdT_C - CpdT_O2_lSc);

    % j) MnO + C -> Mn + CO
    dH_Tj = r_MnO_C * ((dH_CO + dH_MnS - dH_MnO - dH_CS) + ...
        CpdT_Mn + CpdT_CO_lSc - CpdT_MnO - CpdT_C);

    % k) 2MnO + Si -> 2Mn + SiO2
    dH_Tk = r_2MnO_Si * ((dH_SiO2 + dH_SiO2S + 2*dH_MnS - 2*dH_MnO - dH_SiS) + ...
        2*CpdT_Mn + CpdT_SiO2 - CpdT_Si - 2*CpdT_MnO);

    % l) Si + O2 -> SiO2
    dH_Tl = r_Si_O2 * ((dH_SiO2 + dH_SiO2S - dH_SiS) + ...
        CpdT_SiO2 - CpdT_Si - 2*CpdT_O2_lSc);

    % m) 2Cr + 3/2O2 -> Cr2O3
    dH_Tm = r_2Cr_3hO2 * ((dH_Cr2O3 - 2*dH_CrS) + ...
        CpdT_Cr2O3 - 2*CpdT_Cr - 1.5*CpdT_O2_lSc);

    % n) CH4 + 2O2 -> CO2 + 2H2O
%     dH_Tn = -(CH4_inj/M_CH4) * ((dH_CO2 + 2*dH_H2O - dH_CH4) + ...
%         (Cp_CO2 + 2*Cp_H2O - Cp_CH4 - 2*Cp_O2) * (T_gas - 298));

    % o) Graphite to CO2
    dH_To = -(dm_el/M_C) * (dH_CO2 + CpdT_CO2_gas - CpdT_C_gas - CpdT_O2_gas);

    % p) C9H20 + 14O2 -> 9CO2 + 10H2O
    dH_Tp = r_comb * ((9*dH_CO2 + 10*dH_H2O - dH_C9H20) + ...
        9*CpdT_CO2_gas + 10*CpdT_H2O - CpdT_C9H20 - 14*CpdT_O2_gas);
    
    % ======================== Heat Transfer =========================
    
    % --------------- Solid Metal ---------------
    
    % Energy dissipated from the arcs by conduction
    Q_arc = 0.2 * P_arc;

    % Conduction between the solid and liquid metal zones
    Q_lScsSc = min([m_lSc, m_sSc]) * K_therm1 * K_area1 * (T_lSc - T_sSc);

    % Conduction between solid metal and solid slag
    Q_sScsSl = min([m_sSc, m_sSl]) * K_therm2 * K_area2 * (T_sSc - T_sSl);

    % Conduction between solid metal and liquid slag
    Q_sSclSl = min([m_sSc, m_lSl]) * K_therm3 * K_area3 * (T_sSc - T_lSl);
    
    % Convection between solid metal and gas zone
    Q_sScgas = (m_sSc/m_EAF) * K_therm4 * (T_sSc-T_gas) * (1-K_sSclSc);
    
    % Cooling of solid metal to the wall
    Q_sScwater = K_water1 * (T_sSc - T_wall) * (T_sSc/T_melt) * (1 - exp(-(m_sSc/m_EAF)));
    
    % --------------- Liquid Metal ---------------
    
    % Conduction between liquid metal and solid slag
    Q_lScsSl = min([m_lSc, m_sSl]) * K_therm5 * K_area5 * (T_lSc - T_sSl);

    % Conduction between liquid metal and liquid slag
    Q_lSclSl = min([m_lSc, m_lSl]) * K_therm6 * K_area6 * (T_lSc - T_lSl);

    % Convection between the liquid metal and surrounding gas
    Q_lScgas = (m_lSc/m_EAF) * K_therm7 * (T_lSc - T_gas) * K_sSclSc;

    % Cooling of liquid metal to furnace wall
    Q_lScwater = K_water2 * (T_lSc - T_wall) * (T_lSc/T_melt) * (1 - exp(-(m_lSc/m_EAF)));
    
    % --------------- Solid Slag ---------------
    
    % Cooling of solid slag to furnace wall
    Q_sSlwater = K_water3 * (T_sSl - T_wall) * (T_sSl/T_melt) * (1 - exp(-(m_sSl/m_EAF)));
    
    % --------------- Liquid Slag ---------------
    
    % Heat exchange between liquid slag and gas zone
    Q_lSlgas = (m_lSl/m_EAF) * K_therm8 * (T_lSl - T_gas) * K_sSclSc;

    % Energy loss in liquid slag due to cooling
    Q_lSlwater = K_water4 * (T_lSl - T_wall) * (T_lSl/T_melt) * (1 - exp(-(m_lSl/m_EAF)));
    
    % ---------------- Gas Zone -----------------
    
    % Energy received by gas zone from arcs
    Q_arcgas = 0.025 * P_arc;
    
    % Energy loss from gas zone to furnace roof and walls
    Q_gaswater = K_water5*((T_gas - T_roof)*(A1/(A1+A2)) + (T_gas - T_wall)*(A2/(A1+A2)));
    
    % ====================== Reactor Areas ====================
    
    % Areas of roof and wall
    A1 = (pi * r_eafout^2) - (pi * r_hole^2); % roof
    A2 = 2 * pi * r_eafout * h_wall; % wall
    
    % Surface area of sSc and lSc
    A3 = (pi * r_eafout^2) - (pi * (d_coneout/2)^2) + (pi*0.75*d_coneout*sqrt(h_cone + d_coneout/4));
    A4 = pi * r_eafin^2;
    
    % ====================== Slag Foaming =====================
    % All from Fruham 1999
    
    % Slag viscosity
    % Read from the ternary diagram of CaO-FeO-SiO2 system (2.84)
    nu = 0.4;
    
    % Surface tension
    sigma = 0.475; % N/m
    
    % Rate of Oxidation to CO
    r_CO_ox = r_FeO_CL + r_FeO_CD + r_C_hO2 + r_MnO_C;
    
    % Gas Vol. flowrate and Superficial Gas Velocity
    U_g = (r_CO_ox * 8.314 * T_gas) / ((rp + 120000)*A_eaf);
    U_T = (14.55*U_g)/(1-0.089*U_g);
    
    % Bubble Diameter
    % Drag Coeff. Cd assumed unity
    % Bubble velocity 0.7 at U_T = 1.76;
    d_b = (3/(1.14*rho_lSc^2))^(1/3) * ((2*sigma)/(0.7^2));
    
    % Slag Foaming Index
    Xi = 115*nu^1.2/(sigma^0.2*rho_lSl*d_b^0.9); % Taken from Zhang Fruham 1995
    
    % Height change due to foaming
    dh_slag = Xi * U_g;
    
    % Slag Factor
    K_slag = 0.7 * (0.5*tanh(5*(h_lSl+dh_slag)-1.25)+0.5) * (0.5*tanh(3.2*(1-(m_sSc/1000)) - 1.29) + 0.5);
    
    % ======================== View Factor =======================
    
    % =========================
    % ======== Glossary =======
    % =========================
    % ====     1 = Roof    ====
    % ====     2 = Wall    ====
    % ====     3 = sSc     ====
    % ====     4 = lSc     ====
    % ====     5 = arc     ====
    % =========================
    
    % VF_51 Arc -> Roof
    R = r_electrode / r_eafout;
    H1 = (h_electrode + h_arc) / r_eafout;
    H2 = h_electrode / r_eafout;
    a1 = H1^2 + R^2 - 1;
    a2 = H2^2 + R^2 - 1;
    b1 = H1^2 - R^2 + 1;
    b2 = H2^2 - R^2 + 1;
    
    VF_511 = (b1/(8*R*H1)) + (1/(2*pi))*(acos(a1/b1) - (1/(2*H1))*sqrt((a1+2)^2/(R^2)-4) ...
        *acos((a1*R)/b1) - (a1/(2*R*H1))*asin(R));
    VF_512 = (b2/(8*R*H2)) + (1/(2*pi))*(acos(a2/b2) - (1/(2*H2))*sqrt((a2+2)^2/(R^2)-4) ...
        *acos((a2*R)/b2) - (a2/(2*R*H2))*asin(R));
    
    A511 = 2*pi*r_electrode*(h_electrode + h_arc);
    A512 = 2*pi*r_electrode*h_electrode;
    A513 = 2*pi*r_electrode*h_arc;
    
    VF_51 = (1-K_slag)*((VF_511*A511 - VF_512*A512)/A513);
    
    % VF_52 Arc -> Wall
    X = h_cone/r_eafout;
    Y = h_wall/r_eafout;
    L = h_arc/r_eafout;
    R = r_electrode/r_eafout;
    
    aX = X^2 + R^2 - 1;
    bX = X^2 - R^2 + 1;
    FX = (bX/(8*R*X)) + (1/(2*pi))*(acos(aX/bX) - (1/(2*X))*sqrt((aX+2)^2/(R^2)-4) ...
        *acos((aX*R)/bX) - (aX/(2*R*X))*asin(R));
    if isnan(FX)
        FX = 0;
    end
    
    aLX = (L-X)^2 + R^2 - 1;
    bLX = (L-X)^2 - R^2 + 1;
    FLX = (bLX/(8*R*(L-X))) + (1/(2*pi))*(acos(aLX/bLX) - (1/(2*(L-X)))*sqrt((aLX+2)^2/(R^2)-4) ...
        *acos((aLX*R)/bLX) - (aLX/(2*R*(L-X)))*asin(R));
    
    aYXL = (Y+X-L)^2 + R^2 - 1;
    bYXL = (Y+X-L)^2 - R^2 + 1;
    FYXL = (bYXL/(8*R*(Y+X-L))) + (1/(2*pi))*(acos(aYXL/bYXL) - (1/(2*(Y+X-L)))*sqrt((aYXL+2)^2/(R^2)-4) ...
        *acos((aYXL*R)/bYXL) - (aYXL/(2*R*(Y+X-L)))*asin(R));
    
    aXY = (X+Y)^2 + R^2 - 1;
    bXY = (X+Y)^2 - R^2 + 1;
    FXY = (bXY/(8*R*(X+Y))) + (1/(2*pi))*(acos(aXY/bXY) - (1/(2*(X+Y)))*sqrt((aXY+2)^2/(R^2)-4) ...
        *acos((aXY*R)/bXY) - (aXY/(2*R*(X+Y)))*asin(R));
    
    VF_52 = (1-K_slag) * ((X/L)*FX + ((L-X)/L)*(1-FLX) + ((Y+X-L)/L)*FYXL - ((X+Y)/L)*FXY);
    
    % VF_53 Arc -> sSc
    m_charge = m_sSc + m_lSc;
    VF_53 = (1-VF_51-VF_52) * (1-K_sSclSc*(1-(m_sSc/m_charge)));
    
    % VF_54 Arc -> lSc
    VF_54 = (1-VF_51-VF_52) * (K_sSclSc*(1-(m_sSc/m_charge)));
    
    % VF_41 lSc -> Roof
    H = h_wall / r_eafin;
    R2 = r_hole / r_eafin;
    R3 = r_eafout / r_eafin;
    
    VF_41 = (1/2)*(R3^2 - R2^2 - sqrt((1+R3^2+H^2)^2 - 4*R3^2) + ...
        sqrt((1+R2^2+H^2)^2 - 4*R2^2));
    
    % VF_14 Roof -> lSc
    VF_14 = VF_41 * (A4/A1);
    
    % VF_13 Roof -> sSc
%     H = h_wall / r_hole;
%     R2 = r_eafout / r_hole;
%     R3 = d_coneout/2 / r_hole;
%     R4 = d_conein/2 / r_hole;
    
%     VF_131 = 1/(2*(R2^2-1)) * (sqrt((R2^2+R3^2+H^2)^2 - (2*R3*R2)^2) - ...
%         sqrt((R2^2+R4^2+H^2)^2 - (2*R2*R4)^2) + sqrt((1+R4^2+H^2)^2 - (2*R4^2)^2) ...
%         - sqrt((1+R3^2+H^2)^2 - (2*R3^2)^2));
    
    H = h_cone / (d_conein/2);
    R = (d_coneout/2) / (d_conein/2);
    X = 1 + R^2 + H^2;
    
    VF_132 = (2*R^2-X+sqrt(X^2-4*R^2))/(2*sqrt(X-2*R)*(1+R));
    
    VF_13 = VF_132 * VF_41;
    
    % VF_31 sSc -> Roof
    VF_31 = VF_13 * (A1/A3);
    
    
    
    % VF_32 sSc -> Wall
    R = r_eafout / (d_coneout/2);
    H = h_wall / (d_coneout/2);
    
    VF_321 = VF_132 * (1/2) * (1-R^2-H^2+sqrt((1+R^2+H^2)^2-4*R^2));
%     VF_322 = (1/2)*(1 + (1/(R^2-1)) * (H*sqrt(4*R^2+H^2) - ...
%         sqrt((1+R^2+H^2)^2 - 4*R^2)));
    
    VF_32 = VF_321;
    
    % VF_23 Wall -> sSc
    VF_23 = VF_32 * (A3/A2);
    
    % VF_42 lSc -> Wall
    R = r_eafout / r_eafin;
    H = h_wall / r_eafin;
    
    VF_42 = (1/2)*(1-R^2-H^2+sqrt((1+R^2+H^2)^2-4*R^2));
    
    % VF_24 Wall -> lSc
    VF_24 = VF_42 * (A4/A2);
    
    % Inverse
    VF_15 = VF_51 * (A513/A1);
    VF_25 = VF_52 * (A513/A2);
    VF_35 = VF_53 * (A513/A3);
    VF_45 = VF_54 * (A513/A4);
    
    % VF_12 Roof -> Wall
    VF_12 = 1 - VF_13 - VF_14 - VF_15;
    
    % VF_21 Wall -> Roof
    VF_21 = VF_12 * (A1/A2);
    
    % ========================= Radiosity ========================
    Q_arcRAD = 0.75 * P_arc;
    
    % Radiosity of roof
    J_roof = (ep1*sig*T_roof^4/1000 + (1-ep1)*(VF_12*J_wall + VF_13*J_sSc ...
        + VF_14*J_lSc + VF_15*Q_arcRAD));

    % Radiosity of wall
    J_wall = (ep2*sig*T_wall^4/1000 + (1-ep2)*(VF_21*J_wall + VF_23*J_sSc ...
        + VF_24*J_lSc + VF_25*Q_arcRAD));

    % Radiosity of sSc
    J_sSc = (ep3*sig*T_sSc^4/1000 + (1-ep3)*(VF_31*J_roof + VF_32*J_wall ...
        + VF_35*Q_arcRAD));

    % Radiosity of lSc
    J_lSc = (ep4*sig*T_lSc^4/1000 + (1-ep4)*(VF_41*J_roof + VF_42*J_wall ...
        + VF_45*Q_arcRAD));
    
    % ========================= Radiation ========================
    % Radiative heat flow in roof
    Q_roofRAD = A1 * (VF_12*(J_roof-J_wall) + VF_13*(J_roof-J_sSc) ...
        + VF_14*(J_roof-J_lSc)) - VF_51 * Q_arcRAD;

    % Radiative heat flow in wall
    Q_wallRAD = A2 * (VF_21*(J_wall-J_roof) + VF_23*(J_wall-J_sSc) ...
        + VF_24*(J_wall-J_lSc)) - VF_52 * Q_arcRAD;

    % Radiative heat flow in sSc
    Q_sScRAD = A3 * (VF_31*(J_sSc-J_roof) + VF_32*(J_sSc-J_wall)) ...
        - VF_53 * Q_arcRAD;

    % Radiative heat flow in lSc
    Q_lScRAD = A4 * (VF_41*(J_lSc-J_roof) + VF_42*(J_lSc-J_wall)) ...
        - VF_54 * Q_arcRAD;
    
    % ====================== Total Heat Flow =====================
    
    Q_lScchem = (dH_Ta + dH_Tb + dH_Tc + dH_Td + dH_Te + dH_Tf + dH_Tg ...
    + dH_Ti + dH_Tj + dH_Tk + dH_Tl + dH_Tm + dH_Tp);
    
    % Net heat flow in solid steel zone (sSc)
    % CO post combustion and Oxygen burner neglected
    Q_sSc = (Q_arc - dH_Th)*(1-K_sSclSc) + Q_lScsSc ...
    - Q_sScsSl - Q_sSclSl - Q_sScgas - Q_sScwater - Q_sScRAD;
    
    % Net heat flow in liquid metal zone (lSc)
    % CO post combustion and Oxygen burner neglected
    Q_lSc = (Q_arc + dH_Th)*K_sSclSc - Q_lScchem ...
    - Q_lScsSc - Q_lScsSl - Q_lSclSl - Q_lScgas - Q_lScwater - Q_lScRAD;
    
    % Net heat flow in solid slag zone
    Q_sSl = Q_sScsSl + Q_lScsSl - Q_sSlwater;
    
    % Net heat flow in liquid slag zone
    Q_lSl = Q_lSclSl + Q_sSclSl - Q_lSlgas - Q_lSlwater;
    
    % Gas zone energy balance
    Q_gas = Q_arcgas + Q_sScgas + Q_lScgas + Q_lSlgas - Q_gaswater;
    
    
    % ==================== Temperature Change ====================
    
    % Temperature change of sSc
    dT_sSc = (Q_sSc*(1-(T_sSc/T_melt))) / ((m_sSc/M_Fe)*Cp_sSc);

    % Temperature change of lSc
    dT_lSc = (Q_lSc/((m_lSc/M_Fe)*Cp_lSc))/5;
    
    % Temperature change of sSl
    dT_sSl = (Q_sSl*(1-(T_sSl/T_melt))) / ((m_sSl/M_sSl)*Cp_sSl);
    
    % Temperature change of gas
    dT_gas = Q_gas/((m_gas/M_gas)*Cp_gas);
    
    % Temperature change of lSl
    dT_lSl = Q_lSl/((m_lSl/M_lSl)*Cp_lSl);
    if T_lSl < 1820
        dT_lSl = -dT_lSl;
    end

    % Temperature change of roof
    dT_roof = (-Q_roofRAD + (A1/(A1+A2))*Q_gaswater - phi1*Cp_H2O*(T_roof - T_water)) ...
        / (A1*d1*rho*Cp_roof);

    % Temperature change of wall
    dT_wall = (-Q_wallRAD + (A2/(A1+A2))*Q_gaswater - phi2*Cp_H2O*(T_wall - T_water)) ...
        / (A2*d2*rho*Cp_wall);
    
    
    T_lSc = T_lSc + dT_lSc * ts;
    T_lSl = T_lSl + dT_lSl * ts;
    T_gas = T_gas + dT_gas * ts;
    T_roof = T_roof + dT_roof * ts;
    T_wall = T_wall + dT_wall * ts;
    
    % ======================= Phase Change =======================
    
    % Injected Carbon Dissolve Rate
    % Needs to be verified
    dm_CL_melt = (m_CL * T_lSc * Cp_lSc * (T_air/T_melt)) / ...
        (lambda_C + Cp_C * (T_melt - T_air));
    
    % Melt rate of solid metal (kg/s)
    dm_sSc = ((Q_sSc*(T_sSc/T_melt)) / ((lambda_sSc + Cp_sSc*(T_melt - T_sSc))/0.6));
    T_sSc = T_sSc + dT_sSc * ts;
    T_lSc = (T_lSc * m_lSc + T_melt * dm_sSc * ts) / (m_lSc + dm_sSc * ts);
    
    m_Cr_sSc = m_Cr_sSc - (dm_sSc*MX_Cr_sSc*ts);
    m_Cr_lSc = m_Cr_lSc + (dm_sSc*MX_Cr_sSc*ts);
    m_Fe_sSc = m_Fe_sSc - (dm_sSc*MX_Fe_sSc*ts);
    m_Fe_lSc = m_Fe_lSc + (dm_sSc*MX_Fe_sSc*ts);
    m_Mn_sSc = m_Mn_sSc - (dm_sSc*MX_Mn_sSc*ts);
    m_Mn_lSc = m_Mn_lSc + (dm_sSc*MX_Mn_sSc*ts);
    m_P_sSc = m_P_sSc - (dm_sSc*MX_P_sSc*ts);
    m_P_lSc = m_P_lSc + (dm_sSc*MX_P_sSc*ts);
    m_Si_sSc = m_Si_sSc - (dm_sSc*MX_Si_sSc*ts);
    m_Si_lSc = m_Si_lSc + (dm_sSc*MX_Si_sSc*ts);
    m_C_sSc = m_C_sSc - (dm_sSc*MX_C_sSc*ts);
    m_C_lSc = m_C_lSc + (dm_sSc*MX_C_sSc*ts);
    
    % Melt rate of solid slag (kg/s)
    % Melt temperature of 1400C according to https://core.ac.uk/download/pdf/82678298.pdf
    dm_sSl = (Q_sSl*(T_sSl/1673)) / ((lambda_sSl + Cp_sSl*(1673 - T_sSl))/M_sSl);
    T_sSl = T_sSl + dT_sSl * ts;
    T_lSl = (T_lSl * m_lSl + 1673 * dm_sSl * ts) / (m_lSl + dm_sSl * ts);
    
    m_Al2O3_sSl = m_Al2O3_sSl - (dm_sSl*MX_Al2O3_sSl*ts);
    m_Al2O3_lSl = m_Al2O3_lSl + (dm_sSl*MX_Al2O3_sSl*ts);
    m_CaO_sSl = m_CaO_sSl - (dm_sSl*MX_CaO_sSl*ts);
    m_CaO_lSl = m_CaO_lSl + (dm_sSl*MX_CaO_sSl*ts);
    m_MgO_sSl = m_MgO_sSl - (dm_sSl*MX_MgO_sSl*ts);
    m_MgO_lSl = m_MgO_lSl + (dm_sSl*MX_MgO_sSl*ts);
    m_SiO2_sSl = m_SiO2_sSl - (dm_sSl*MX_SiO2_sSl*ts);
    m_SiO2_lSl = m_SiO2_lSl + (dm_sSl*MX_SiO2_sSl*ts);
    
    % ===================== Geometry Change ======================
    
    % height of liquid metal
    h_lSc = (m_lSc / rho_lSc) / A_bath;

    % height of liquid slag
    h_lSl = (m_lSl / rho_lSl) / A_bath;
    
    h_sSc2 = 0;
    d_conein = r_eafin*2;
    d_coneout = r_eafout*2;

    V_sSc = m_sSc / rho_sSc;
    h_cone = V_sSc / (pi*r_eafout^2-(1/3)*pi*(r_eafin^2 + r_eafin*r_eafout + r_eafout^2));
    h_sSc1 = h_cone;

    % Arc height
    h_arc = h_eafup - h_electrode - (h_lSc + h_lSl);

    % height of wall
    h_wall = h_eafup + h_eaflow - h_sSc1 - h_lSc - h_lSl;
    
    % Exposure Coeff.
    K_sSclSc = 0.5*tanh(5*(h_lSc-h_sSc1-h_sSc2+h_cone)) + 0.5;
    
    % =================== Reaction Mass Change ===================
    
    % Injected carbon dissolution
    m_CL = m_CL - (dm_CL_melt * ts);
    m_C_lSc = m_C_lSc + (dm_CL_melt * ts);
    
    % ---------- Decarburization ---------
    % FeO + C -> Fe + CO
    
    % By injected carbon
    m_FeO_lSl = m_FeO_lSl - (r_FeO_CL * M_FeO * ts);
    m_CL = m_CL - (r_FeO_CL * M_C * ts);
    m_Fe_lSc = m_Fe_lSc + (r_FeO_CL * M_Fe * ts);
    m_CO = m_CO + (r_FeO_CL * M_CO * ts);
    
    % By dissolved carbon
    m_FeO_lSl = m_FeO_lSl - (r_FeO_CD * M_FeO * ts);
    m_C_lSc = m_C_lSc - (r_FeO_CD * M_C * ts);
    m_Fe_lSc = m_Fe_lSc + (r_FeO_CD * M_Fe * ts);
    m_CO = m_CO + (r_FeO_CD * M_CO * ts);
    
    % ---------- Carbon Oxidation -----------
    
    % To carbon monoxide
    % C + 1/2 O2 -> CO
    m_C_lSc = m_C_lSc - (r_C_hO2 * M_C * ts);
    m_CO = m_CO + (r_C_hO2 * M_CO * ts);
    m_O2 = m_O2 - 0.5*(r_C_hO2 * M_O2 * ts);
    
    % To carbon dioxide
    % C + O2 -> CO2
    m_C_lSc = m_C_lSc - (r_C_O2 * M_C * ts);
    m_CO2 = m_CO2 + (r_C_O2 * M_CO2 * ts);
    m_O2 = m_O2 - (r_C_O2 * M_O2 * ts);
    
    % -------- MnO Decarburization ----------
    
    % MnO + C -> Mn + CO
    
    m_MnO_lSl = m_MnO_lSl - (r_MnO_C * M_MnO * ts);
    m_C_lSc = m_C_lSc - (r_MnO_C * M_C * ts);
    m_Mn_lSc = m_Mn_lSc + (r_MnO_C * M_Mn * ts);
    m_CO = m_CO + (r_MnO_C * M_MnO * ts);
    
    % ----------- Desiliconization ----------
    
    % 2FeO + Si -> 2Fe + SiO2
    
    m_FeO_lSl = m_FeO_lSl - 2*(r_2FeO_Si * M_FeO * ts);
    m_Si_lSc = m_Si_lSc - (r_2FeO_Si * M_Si * ts);
    m_Fe_lSc = m_Fe_lSc + 2*(r_2FeO_Si * M_Fe * ts);
    m_SiO2_lSl = m_SiO2_lSl + (r_2FeO_Si * M_SiO2 * ts);
    
    % --------- Silicon Oxidation ---------
    
    % Si + O2 -> SiO2
    
    m_Si_lSc = m_Si_lSc - (r_Si_O2 * M_Si * ts);
    m_SiO2_lSl = m_SiO2_lSl + (r_Si_O2 * M_SiO2 * ts);
    m_O2 = m_O2 - (r_Si_O2 * M_O2 * ts);
    
    % ------- Si reaction with MnO --------
    
    % 2MnO + Si -> 2Mn + SiO2
    
    m_MnO_lSl = m_MnO_lSl - 2*(r_2MnO_Si * M_MnO * ts);
    m_Si_lSc = m_Si_lSc - (r_2MnO_Si * M_Si * ts);
    m_Mn_lSc = m_Mn_lSc + 2*(r_2MnO_Si * M_Mn * ts);
    m_SiO2_lSl = m_SiO2_lSl + (r_2MnO_Si * M_SiO2 * ts);
    
    % ------- Mn reaction with FeO ---------
    
    % Mn + FeO -> MnO + Fe
    
    m_Mn_lSc = m_Mn_lSc - (r_FeO_Mn * M_Mn * ts);
    m_FeO_lSl = m_FeO_lSl - (r_FeO_Mn * M_FeO * ts);
    m_MnO_lSl = m_MnO_lSl + (r_FeO_Mn * M_MnO * ts);
    m_Fe_lSc = m_Fe_lSc + (r_FeO_Mn * M_Fe * ts);
    
    % ------- Cr reaction with FeO ---------
    
    % 3FeO + 2Cr -> 3Fe + Cr2O3
    
    m_FeO_lSl = m_FeO_lSl - 3*(r_3FeO_2Cr * M_FeO * ts);
    m_Cr_lSc = m_Cr_lSc - 2*(r_3FeO_2Cr * M_Cr * ts);
    m_Fe_lSc = m_Fe_lSc + 3*(r_3FeO_2Cr * M_Fe * ts);
    m_Cr2O3_lSl = m_Cr2O3_lSl + (r_3FeO_2Cr * M_Cr2O3 * ts);
    
    % --------- Chromium Oxidation ---------
    
    % 2Cr + 3/2O2 -> Cr2O3
    
    m_Cr_lSc = m_Cr_lSc - 2*(r_2Cr_3hO2 * M_Cr * ts);
    m_Cr2O3_lSl = m_Cr2O3_lSl + (r_2Cr_3hO2 * M_Cr2O3 * ts);
    m_O2 = m_O2 - 1.5*(r_2Cr_3hO2 * M_O2 * ts);
    
    % -------- Phosphorus Oxidation --------
    
    % 5FeO + 2P -> 5Fe + P2O5
    
    m_FeO_lSl = m_FeO_lSl - 5*(r_5FeO_2P * M_FeO * ts);
    m_P_lSc = m_P_lSc - 2*(r_5FeO_2P * M_P * ts);
    m_Fe_lSc = m_Fe_lSc + 5*(r_5FeO_2P * M_Fe * ts);
    m_P2O5_lSl = m_P2O5_lSl + (r_5FeO_2P * M_P2O5 * ts);
    
    % ------------- Fe Oxidation -------------
    
    % Fe + 1/2O2 -> FeO
    
    m_Fe_lSc = m_Fe_lSc - (r_Fe_hO2 * M_Fe * ts);
    m_FeO_lSl = m_FeO_lSl + (r_Fe_hO2 * M_FeO * ts);
    
    % ------------- Combustion ---------------
    
    % C9H20 + 14O2 -> 9CO2 + 10H2O
    
    m_comb_sSc = m_comb_sSc - (r_comb * M_C9H20 * ts);
    m_O2 = m_O2 - 14*(r_comb * M_O2 * ts);
    m_CO = m_CO + 9*(r_comb * M_CO2 * ts);
    m_H2O = m_H2O + 10*(r_comb * M_H2O * ts);
    
    % ----------- Pose Combustion ------------
    
    % CO + 1/2O2 -> CO2
    
    m_O2 = m_O2 + (O2_post - r_post*M_O2) * ts;
    m_CO = m_CO - 2*(r_post * M_CO * ts);
    m_CO2 = m_CO2 + (r_post * M_CO2 * ts);
    
    % --------- Electrode Oxidation ----------
    
    % C + O2 -> CO2
    
    m_O2 = m_O2 - ((dm_el*M_O2)/M_C) * ts;
    m_CO2 = m_CO2 + ((dm_el*M_CO2)/M_C) * ts;
    
    % ======================= Material Addition ======================
    
    % ------------ Solid Slag ------------
    m_CaO_sSl = m_CaO_sSl + (slg_add*ts) * MX_CaO_slg;
    m_MgO_sSl = m_MgO_sSl + (slg_add*ts) * MX_MgO_slg;
    m_SiO2_sSl = m_SiO2_sSl + (slg_add*ts) * MX_SiO2_slg;
    m_Al2O3_sSl = m_Al2O3_sSl + (slg_add*ts) * MX_Al2O3_slg;
    
    T_sSl = (T_sSl * m_sSl + T_slg * slg_add * ts) / (m_sSl + slg_add * ts);
    
    % ------ Extra Material Addition -----
    
    % Oxygen Lance
    m_O2 = m_O2 + O2_lance * ts;
    
    T_gas = (T_gas * m_gas + T_air * O2_lance * ts) / (m_gas + O2_lance * ts);
    
    % Carbon injection
    m_CL = m_CL + (C_inj * ts);
    
    % Manganese Injection
    m_Mn_sSc = m_Mn_sSc + (Mn_inj * ts);
    
    % ======================== Take out ========================
%     thres_gas = (121590*V_gas) / (R * T_gas);
% 
%     if XM_gas > thres_gas
%         gas_out = XM_gas - thres_gas;
%         H2O_out = gas_out * X_H2O;
%         CO2_out = gas_out * X_CO2;
%         CO_out = gas_out * X_CO;
%         O2_out = gas_out * X_O2;
%         
%         m_H2O = m_H2O - (H2O_out*M_H2O);
%         m_CO2 = m_CO2 - (CO2_out*M_CO2);
%         m_CO = m_CO - (CO_out*M_CO);
%         m_O2 = m_O2 - (O2_out*M_O2);
%         
%         T_gas = (T_gas * XM_gas + T_air * gas_out) / (XM_gas + gas_out);
%     end

    % Off gas venting
    m_CO = m_CO - (hd*u1*MX_CO)/(k_U*u2+hd) * ts;
    m_CO2 = m_CO2 - (hd*u1*MX_CO2)/(k_U*u2+hd) * ts;
    m_O2 = m_O2 - (hd*u1*MX_O2)/(k_U*u2+hd) * ts;
    m_H2O = m_H2O - (hd*u1*MX_H2O)/(k_U*u2+hd) * ts;
    
    dm_CO = r_FeO_CL * M_CO + r_FeO_CD * M_CO + r_C_hO2 * M_CO + r_MnO_C * M_MnO ...
        + 9*r_comb * M_CO2 - 2*r_post * M_CO;
    dm_CO2 = r_C_O2 * M_CO2 + 9*r_comb * M_CO2 + r_post * M_CO2;
    dm_O2 = O2_lance -0.5*r_C_hO2 * M_O2 - r_C_O2 * M_O2 - r_Si_O2 * M_O2 ...
        -1.5*r_2Cr_3hO2 * M_O2 - 14*r_comb * M_O2 + (O2_post - r_post*M_O2);
    dm_H2O = 10*r_comb * M_H2O;
    
    rp = (R*T_gas/V_gas)*(dm_CO/M_CO + dm_CO2/M_CO2 + dm_O2/M_O2 + dm_H2O/M_H2O) ...
        + (R*dT_gas/V_gas)*(m_CO/M_CO + m_CO2/M_CO2 + m_O2/M_O2 + m_H2O/M_H2O);
    
    % ===================== For Graph =======================
    if mod(step, 200) == 0
        gas_temp(floor(step*ts*5)) = T_gas;
        sSc_temp(floor(step*ts*5)) = T_sSc;
        sSl_temp(floor(step*ts*5)) = T_sSl;
        lSc_temp(floor(step*ts*5)) = T_lSc;
        lSl_temp(floor(step*ts*5)) = T_lSl;
%         steel_Fe(step*ts*5) = MX_Fe_lSc;
        m_solid(floor(step*ts*5)) = m_sSc;
        m_liquid(floor(step*ts*5)) = m_lSc;
        m_solid_slag(floor(step*ts*5)) = m_sSl;
        m_liquid_slag(floor(step*ts*5)) = m_lSl;
%         rel_pres(step*ts*5) = rp;
    end
    
end

% Graph generation
time = linspace(1, secs*5, secs*5);

figure
plot(time, gas_temp)
hold on
plot(time, sSc_temp)
plot(time, sSl_temp)
plot(time, lSc_temp)
plot(time, lSl_temp)

legend('Gas', 'Solid Metal', 'Solid Slag', 'Liquid Metal', 'Liquid Slag')
hold off

% figure
% plot(time, steel_Fe)

% figure
% plot(time, rel_pres)

figure
plot(time, m_solid)
hold on
plot(time, m_solid_slag)
plot(time, m_liquid)
plot(time, m_liquid_slag)
legend('solid', 'solid slag', 'liquid', 'liquid slag')
hold off