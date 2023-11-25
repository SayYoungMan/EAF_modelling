clear
clc

% ========================== Control Panel ===========================

% ---------- Time Settings ----------

% Time slice
ts = 1/1000; % 10^-3 s

% Total operating time in s
secs = 18000;

% Takeout interval in s
out = 900;

% ---------- DRI Settings -----------

% Temperature of DRI in K
T_DRI = 1213.2;

% DRI mass addition rate in kg/s
DRI_add = 344912/3600;

% DRI mass fraction
MX_Fe_DRI = 0.89659;
MX_C_DRI = 0.00588;
MX_SiO2_DRI = 0.05772;
MX_Al2O3_DRI = 0.03665;
MX_CaO_DRI = 0.00137;
MX_MgO_DRI = 0.0008;
MX_MnO_DRI = 0.0001;
MX_P2O5_DRI = 0.00188;

% --------- Scrap Settings -----------

% Temperature of scrap in K
T_scr = 300;

% DRI mass addition rate in kg/s
scr_add = 80000/3600;

% DRI mass fraction
MX_Fe_scr = 0.9705;
MX_C_scr = 0.004;
MX_Si_scr = 0.006;
MX_Cr_scr = 0.002;
MX_P_scr = 0.0005;
MX_Mn_scr = 0.006;
MX_comb_scr = 0.011;

% ---------- Slag Settings -----------

% Temperature of slag in K
T_slg = 300;

% Slag mass addition rate in kg/s
slg_add = 1.5;

% Slag mass fraction
MX_CaO_slg = 0.573;
MX_MgO_slg = 0.415;
MX_SiO2_slg = 0.007;
MX_Al2O3_slg = 0.005;

% --------- Reactor Geometry ---------
r_eafout = 3.72;
r_eafin = 3.62;
r_hole = 0.3;
r_electrode = 0.3;
h_eafup = 2.9;
h_eaflow = 1.4;
h_electrode = 1.0;
d1 = 0.30;
d2 = 0.45;

% ---------- Other Settings ----------

% Carbon Injection Rate (kg/s)
C_inj = 0.95;

% Ferro-Manganese Injection Rate (kg/s)
FM_inj = 1.25;

MX_Mn_FM = 0.78;
MX_C_FM = 0.07;
MX_P_FM = 0.002;
MX_Si_FM = 0.003;
MX_Fe_FM = 0.145;

% Oxygen Lance Rate (kg/s)
O2_lance = 2.8;

% O2 for post combustion (kg/s)
O2_post = 0.8;

% Power of arc (kW)
P_arc = 30000;

% EAF mass capacity (kg)
m_EAF = 220000;

% Cooling water flowrate (mol/s)
phi1 = 80/0.018;
phi2 = 150/0.018;

% ======================= Initial Parameters =========================

% ------------- Initial mass (kg) --------------

% Solid metal initial mass
m_Fe_sSc = 1639;
m_C_sSc = 11.5;
m_Cr_sSc = 0.42;
m_Mn_sSc = 19.84;
m_P_sSc = 0.14;
m_SiO2_sSc = 84.2;
m_Al2O3_sSc = 53.5;
m_CaO_sSc = 1.67;
m_MgO_sSc = 0.98;
m_MnO_sSc = 0.123;
m_P2O5_sSc = 2.28;
m_Si_sSc = 1.7;
m_comb_sSc = 3.48;

m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_SiO2_sSc + ...
    m_Al2O3_sSc + m_CaO_sSc + m_MgO_sSc + m_MnO_sSc + m_P2O5_sSc + m_Si_sSc + m_comb_sSc;

% Liquid metal initial mass
m_Fe_lSc = 64595;
m_C_lSc = 321;
m_Cr_lSc = 26.5;
m_Mn_lSc  = 502;
m_P_lSc = 4.5;
m_Si_lSc = 134.6;

m_lSc = m_Fe_lSc + m_C_lSc + m_Cr_lSc + m_Mn_lSc + m_P_lSc + m_Si_lSc;

% Solid slag initial mass
m_CaO_sSl = 56.84;
m_MgO_sSl = 41.2;
m_SiO2_sSl = 2;
m_Al2O3_sSl = 1.4;

% Liquid slag initial mass
m_SiO2_lSl = 3229;
m_Al2O3_lSl = 2153;
m_CaO_lSl = 1119;
m_MgO_lSl = 795;
m_MnO_lSl = 960;
m_P2O5_lSl = 115.3;
m_Cr2O3_lSl = 2.19;
m_FeO_lSl = 3195;

m_lSl = m_SiO2_lSl + m_Al2O3_lSl + m_CaO_lSl + m_MgO_lSl + m_MnO_lSl + ...
    m_P2O5_lSl + m_Cr2O3_lSl + m_FeO_lSl;

% Gas initial mass
m_H2O = 182;
m_O2 = 414;
m_CO = 1388;
m_CO2 = 190;

% Initial mass of injected carbon
m_CL = 0.792;

% ------------- Initial temp. (K) --------------

T_sSc = 1270;
T_lSc = 2001;
T_sSl = 1051;
T_lSl = 1939;
T_gas = 1755;
T_wall = 312;
T_roof = 322;

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
M_MgO = 0.040304;
M_C9H20 = 0.1282;
M_gas = 0.035;
M_Al2O3 = 0.10196;
M_Al = 0.02698;
M_H2O = 0.01802;
M_sSl = 0.0484;
M_lSl = 0.0509;

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
K_therm5 = 0.2;
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
u1 = 20;
u2 = 0.3;

V_gas = 45;

% Electrode
R_tip = 0.02;
R_side = 10;
A_side = 35;
I_arc = 30;

% ------------------------- Arrays for graph ------------------------
gas_temp = zeros(1,secs);
sSc_temp = zeros(1,secs);
sSl_temp = zeros(1,secs);
lSc_temp = zeros(1,secs);
lSl_temp = zeros(1,secs);
steel_Fe = zeros(1,secs);
m_solid = zeros(1,secs);
m_liquid = zeros(1,secs);
m_solid_slag = zeros(1,secs);
m_liquid_slag = zeros(1,secs);
rel_pres = zeros(1,secs);

for step = 1:secs/ts
    
    % ====================== Mole, Mass Fraction =====================
    
    % ------------ Solid Metal ------------
    
    % Total mass of solid metal [kg]
    m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_SiO2_sSc ...
        + m_Al2O3_sSc + m_CaO_sSc + m_MgO_sSc + m_MnO_sSc + m_Si_sSc + ...
        m_comb_sSc + m_P2O5_sSc;
    
    % Total mole of solid metal [mol]
    XM_sSc = (m_Fe_sSc/M_Fe) + (m_C_sSc/M_C) + (m_Cr_sSc/M_Cr) + (m_Mn_sSc/M_Mn) ... 
        + (m_P_sSc/M_P) + (m_SiO2_sSc/M_SiO2) + (m_Al2O3_sSc/M_Al2O3) + ...
        (m_CaO_sSc/M_CaO) + (m_MgO_sSc/M_MgO) + (m_MnO_sSc/M_MnO) + ...
        (m_Si_sSc/M_Si) + (m_comb_sSc/M_C9H20) + (m_P2O5_sSc/M_P2O5);
    
    % Mole fractions of compounds in solid metal
    X_Fe_sSc = (m_Fe_sSc/M_Fe) / XM_sSc;
    X_C_sSc = (m_C_sSc/M_C) / XM_sSc;
    X_Cr_sSc = (m_Cr_sSc/M_Cr) / XM_sSc;
    X_Mn_sSc = (m_Mn_sSc/M_Mn) / XM_sSc;
    X_P_sSc = (m_P_sSc/M_P) / XM_sSc;
    X_SiO2_sSc = (m_SiO2_sSc/M_SiO2) / XM_sSc;
    X_Al2O3_sSc = (m_Al2O3_sSc/M_Al2O3) / XM_sSc;
    X_CaO_sSc = (m_CaO_sSc/M_CaO) / XM_sSc;
    X_MgO_sSc = (m_MgO_sSc/M_MgO) / XM_sSc;
    X_MnO_sSc = (m_MnO_sSc/M_MnO) / XM_sSc;
    X_Si_sSc = (m_Si_sSc/M_Si) / XM_sSc;
    X_comb_sSc = (m_comb_sSc/M_C9H20) / XM_sSc;
    X_P2O5_sSc = (m_P2O5_sSc/M_P2O5) / XM_sSc;
    
    % Mass fractions of compounds in solid metal
    MX_Fe_sSc = m_Fe_sSc / m_sSc;
    MX_C_sSc = m_C_sSc / m_sSc;
    MX_Cr_sSc = m_Cr_sSc / m_sSc;
    MX_Mn_sSc = m_Mn_sSc / m_sSc;
    MX_P_sSc = m_P_sSc / m_sSc;
    MX_SiO2_sSc = m_SiO2_sSc / m_sSc;
    MX_Al2O3_sSc = m_Al2O3_sSc / m_sSc;
    MX_CaO_sSc = m_CaO_sSc / m_sSc;
    MX_MgO_sSc = m_MgO_sSc / m_sSc;
    MX_MnO_sSc = m_MnO_sSc / m_sSc;
    MX_Si_sSc = m_Si_sSc / m_sSc;
    MX_comb_sSc = m_comb_sSc / m_sSc;
    MX_P2O5_sSc = m_P2O5_sSc / m_sSc;
    
    % ----------- Liquid Metal ------------
    
    % Total mass of liquid metal [kg]
    m_lSc = m_Fe_lSc + m_C_lSc + m_Cr_lSc + m_Mn_lSc + m_P_lSc + m_Si_lSc;
    
    % Total mole of liquid metal [mol]
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
    
    % Total mass of solid slag [kg]
    m_sSl = m_CaO_sSl + m_MgO_sSl + m_SiO2_sSl + m_Al2O3_sSl;
    
    % Total mole of solid slag [mol]
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
    
    % Total mass of liquid slag [kg]
    m_lSl = m_SiO2_lSl + m_Al2O3_lSl + m_CaO_lSl + m_MgO_lSl + m_MnO_lSl ...
        + m_P2O5_lSl + m_Cr2O3_lSl + m_FeO_lSl;
    
    % Total mole of liquid slag [mol]
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
    
    % Molar mass of liquid slag [kg/mol]
    M_lSl = X_SiO2_lSl*M_SiO2 + X_Al2O3_lSl*M_Al2O3 + X_CaO_lSl*M_CaO + X_MgO_lSl*M_MgO + ...
        X_MnO_lSl*M_MnO + X_P2O5_lSl*M_P2O5 + X_Cr2O3_lSl*M_Cr2O3 + X_FeO_lSl*M_FeO;
    
    % --------------- Gas -----------------
    
    % Total mass of gas
    m_gas = m_H2O + m_CO + m_CO2 + m_O2;
    
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
    
    % Activity calculation
    gamma_FeO = 10 ^ (1262/T_lSl - 1.1302*X_FeO_lSl + 0.96*X_SiO2_lSl + ...
        0.123*X_CaO_lSl - 0.4198); % Basu et al. 2008
    
    % Equilibrium constant
    K_FC = 10 ^ (-5730/T_lSl + 5.096); % Turkdogan 1996
    
    % Partial pressure
    p_CO = p_gas * X_CO;
    
    % Equilibrium carbon mass/ mole fraction
    k_C = p_CO / (K_FC * (gamma_FeO/(M_FeO*1.65)));
    k_XC = k_C * ((M_lSl * M_Fe) / (M_FeO*M_C*100^2));
    Xeq_C = k_XC * ((m_lSl*M_FeO)/(m_FeO_lSl*M_lSl) + (m_SiO2_lSl*M_FeO) ...
        /(m_FeO_lSl*M_SiO2) + 1); % Bekkar 1999
    if isnan(Xeq_C)
        Xeq_C = 0;
    end
    % Rate of decarburization of dissolved C
    r_FeO_CD = (kd_CD * (X_C_lSc - Xeq_C)) / M_C; % Logar 2012
    
    % ---------- Carbon Oxidation -----------
    
    % To carbon monoxide
    % C + 1/2 O2 -> CO
    r_C_hO2 = (kd_C1 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO) / M_C; % Logar 2012
    if r_C_hO2 > (O2_lance * K_O2CO) / M_C
        r_C_hO2 = (O2_lance * K_O2CO) / M_C;
    end
    
    % To carbon dioxide
    % C + O2 -> CO2
    r_C_O2 = (kd_C2 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO2) / M_C; % Logar 2012
    if r_C_O2 > (O2_lance * K_O2CO2) / M_C
        r_C_O2 = (O2_lance * K_O2CO2) / M_C;
    end
    
    % -------- MnO decarburization ----------
    
    % MnO + C -> Mn + CO
    
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
    
    if isnan(Xeq_Si)
        Xeq_Si = 0;
    end

    % Rate of reaction
    r_2FeO_Si = (kd_Si1 * (X_Si_lSc - Xeq_Si)) / M_Si;
    
    % --------- Silicon Oxidation ---------
    
    % Si + O2 -> SiO2
    
    r_Si_O2 = (kd_Si2 * (X_Si_lSc - Xeq_Si) * O2_lance * K_O2SiO2) / M_Si;
    if r_Si_O2 > (O2_lance * K_O2SiO2) / M_Si
        r_Si_O2 = (O2_lance * K_O2SiO2) / M_Si;
    end
    
    
    % ------- Si reaction with MnO --------
    
    % 2MnO + Si -> 2Mn + SiO2
    
    % Reaction equilibrium constant
    kX_Mn2 = 10^(2.8*B3-1.16) * ((M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2)); % Logar 2012
    
    % Equilibrium MnO mole fraction
    Xeq_MnO2 = sqrt((X_Mn_lSc^2 * X_SiO2_lSl) / (X_Si_lSc * kX_Mn2)); % Logar 2012
    
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
    r_3FeO_2Cr = (2 * kd_Cr1 * (X_Cr_lSc - Xeq_Cr)) / M_Cr;
    
    % --------- Chromium Oxidation ---------
    
    % 2Cr + 3/2O2 -> Cr2O3
    
    r_2Cr_3hO2 = (2 * kd_Cr2 * (X_Cr_lSc - Xeq_Cr) * O2_lance * K_O2Cr2O3) / M_Cr;
    if r_2Cr_3hO2 > (O2_lance * K_O2Cr2O3) / M_Cr
        r_2Cr_3hO2 = (O2_lance * K_O2Cr2O3) / M_Cr;
    end
    
    % ----- Phosphorus reaction with FeO -----
    
    % 5FeO + 2P -> 5Fe + P2O5
    
    partition = 10 ^ (1.97*X_CaO_lSl + 2.0*X_FeO_lSl - 2.04*X_SiO2_lSl + 6713/T_lSl - 1.84); % Basu, 2007
    eq_P = (2 * (m_P2O5_lSl/M_P2O5) * M_P) / partition;
    Xeq_P = eq_P / XM_lSc;
    
    if isnan(Xeq_P)
        Xeq_P = 0;
    end
    r_5FeO_2P = (2 * kd_P * (X_P_lSc - Xeq_P)) / M_P;
    
    % ------------- Fe Oxidation -------------
    
    % Fe + 1/2O2 -> FeO
    
    r_Fe_hO2 = (O2_lance * K_O2FeO) / M_Fe;
    
    % ------------- Combustion ---------------
    
    % C9H20 + 14O2 -> 9CO2 + 10H2O
    
    r_comb = (kd_comb * m_comb_sSc * (T_sSc/T_melt)) / M_C9H20;
    
    % No reaction if no O2 avilable
    if r_comb > (m_O2/M_O2)
        r_comb = 0;
    end
    
    % ----------- Post Combustion ------------
    
    % CO + 1/2O2 -> CO2
    
    K_mCO = 0.9;
    
    r_post = (O2_post * K_mCO) / M_O2;
    
    % --------- Electrode Oxidation ----------
    
    % C + O2 -> CO2
    
    dm_el = 3*((R_tip * (I_arc^2/3600)) + (R_side * (A_side/3600)));
    
    % No reaction if no O2 avilable
    if dm_el > m_O2
        dm_el = 0;
    end
    
    % =================== Reaction Heat Transfer ====================
    
    % ------------ Dynamic heat capacity calculation ----------
    % Units in kJ/mol K
    
    if T_gas < 500
        t_H2O = (T_gas-298)/1000;
        CpdT_H2O = -203.6060*t_H2O + (1523.290/2)*(t_H2O^2) + (-3196.413/3)*(t_H2O^3) ...
        + (2474.455/4)*(t_H2O^4) - (3.855326)/t_H2O;
    elseif (500 <= T_gas) && (T_gas < 1700)
        t_H2O = (500-298)/1000;
        CpdT_H2O1 = -203.6060*t_H2O + (1523.290/2)*(t_H2O^2) + (-3196.413/3)*(t_H2O^3) ...
        + (2474.455/4)*(t_H2O^4) - (3.855326)/t_H2O;
        t_H2O_2 = (T_gas-298)/1000;
        CpdT_H2O2 = 30.0920*t_H2O_2 + (6.832514/2)*(t_H2O_2^2) + (6.793435/3)*(t_H2O_2^3) ... 
        + (-2.534480/4)*(t_H2O_2^4) - (0.082139)/t_H2O_2;
        CpdT_H2O = CpdT_H2O1 + CpdT_H2O2; 
    else
        t_H2O = (500-298)/1000;
%         CpdT_H2O1 = -203.6060*t_H2O + (1523.290/2)*(t_H2O^2) + (-3196.413/3)*(t_H2O^3) ...
%         + (2474.455/4)*(t_H2O^4) - (3.855326)/t_H2O;
        CpdT_H2O1 = 79 * t_H2O;
        t_H2O_2 = (1700-500)/1000;
        CpdT_H2O2 = 30.0920*t_H2O_2 + (6.832514/2)*(t_H2O_2^2) + (6.793435/3)*(t_H2O_2^3) ... 
        + (-2.534480/4)*(t_H2O_2^4) - (0.082139)/t_H2O_2;
        t_H2O_3 = (T_gas-1700)/1000;
%         CpdT_H2O3 = 41.96426*t_H2O_3 + (8.622053/2)*(t_H2O_3^2) + (-1.499780/3)*(t_H2O_3^3) ...
%         + (0.098119/4)*(t_H2O_3^4) - (-11.15764)/t_H2O_3;
        CpdT_H2O3 = 49.75 * t_H2O_3;
        CpdT_H2O = CpdT_H2O1 + CpdT_H2O2 + CpdT_H2O3; 
    end
    
    t_C = (T_lSc-298)/1000;
    CpdT_C = 21.17510*t_C + (-0.812428/2)*(t_C^2) + (0.448537/3)*(t_C^3) + (-0.043256/4)*(t_C^4) ...
        - (-0.013103)/t_C;
        
    t_C_gas = (T_gas-298)/1000;
    CpdT_C_gas = 21.17510*t_C_gas + (-0.812428/2)*(t_C_gas^2) + (0.448537/3)*(t_C_gas^3) ...
    + (-0.043256/4)*(t_C_gas^4) - (-0.013103)/t_C_gas;
    
    if T_lSc < 1650
        t_FeO = (T_lSc-298)/1000;
        CpdT_FeO = 45.75120*t_FeO + (18.78553/2)*(t_FeO^2) + (-5.952201/3)*(t_FeO^3) ...
        + (0.852779/4)*(t_FeO^4) - (-0.081265)/t_FeO;
    else
        t_FeO = (1650-298)/1000;
        CpdT_FeO1 = 45.75120*t_FeO + (18.78553/2)*(t_FeO^2) + (-5.952201/3)*(t_FeO^3) ...
        + (0.852779/4)*(t_FeO^4) - (-0.081265)/t_FeO;
        t_FeO_2 = (T_lSc-1650)/1000;
        CpdT_FeO2 = 68.1992*t_FeO_2;
        CpdT_FeO = CpdT_FeO1 + CpdT_FeO2;
    end
    
    t_Fe = (T_lSc-298)/1000;
    CpdT_Fe = 23.97449*t_Fe + (8.367750/2)*(t_Fe^2) + (0.000277/3)*(t_Fe^3) + (-0.000086/4)*(t_Fe^4) ...
            - (-0.000005)/t_Fe;
    
    if T_gas < 700
        t_O2_gas = (T_gas-298)/1000;
        CpdT_O2_gas = 31.32234*t_O2_gas + (-20.23532/2)*(t_O2_gas^2) + (57.86644/3)*(t_O2_gas^3) ...
        + (-36.50624/4)*(t_O2_gas^4) - (-0.007374)/t_O2_gas;
    else
        t_O2_gas = (700-298)/1000;
        CpdT_O2_gas1 = 31.32234*t_O2_gas + (-20.23532/2)*(t_O2_gas^2) + (57.86644/3)*(t_O2_gas^3) ...
        + (-36.50624/4)*(t_O2_gas^4) - (-0.007374)/t_O2_gas;
        t_O2_gas2 = (T_gas-700)/1000;
        CpdT_O2_gas2 = 30.03235*t_O2_gas2 + (8.772972/2)*(t_O2_gas2^2) + (-3.988133/3)*(t_O2_gas2^3) ...
        + (0.788313/4)*(t_O2_gas2^4) - (-0.741599)/t_O2_gas2;
        CpdT_O2_gas = CpdT_O2_gas1 + CpdT_O2_gas2;
    end
    
    if T_lSc < 700
        t_O2_lSc = (T_lSc-298)/1000;
        CpdT_O2_lSc = 31.32234*t_O2_lSc + (-20.23532/2)*(t_O2_lSc^2) + (57.86644/3)*(t_O2_lSc^3) ...
        + (-36.50624/4)*(t_O2_lSc^4) - (-0.007374)/t_O2_lSc;
    else
        t_O2_lSc1 = (700-298)/1000;
        CpdT_O2_lSc1 = 31.32234*t_O2_lSc1 + (-20.23532/2)*(t_O2_lSc1^2) + (57.86644/3)*(t_O2_lSc1^3) ...
        + (-36.50624/4)*(t_O2_lSc1^4) - (-0.007374)/t_O2_lSc1;
        t_O2_lSc2 = (T_lSc-700)/1000;
        CpdT_O2_lSc2 = 30.03235*t_O2_lSc2 + (8.772972/2)*(t_O2_lSc2^2) + (-3.988133/3)*(t_O2_lSc2^3) ...
        + (0.788313/4)*(t_O2_lSc2^4) - (-0.741599)/t_O2_lSc2;
        CpdT_O2_lSc = CpdT_O2_lSc1 + CpdT_O2_lSc2;
    end
    
    if T_gas < 1300
        t_CO_gas = (T_gas-298)/1000;
        CpdT_CO_gas = 25.56759*t_CO_gas + (6.096130/2)*(t_CO_gas^2) + (4.054656/3)*(t_CO_gas^3) ...
        + (-2.671301/4)*(t_CO_gas^4) - (0.131021)/t_CO_gas;
    else
        t_CO_gas = (1300-298)/1000;
        CpdT_CO_gas_1 = 25.56759*t_CO_gas + (6.096130/2)*(t_CO_gas^2) + (4.054656/3)*(t_CO_gas^3) ...
        + (-2.671301/4)*(t_CO_gas^4) - (0.131021)/t_CO_gas;
        t_CO_gas_2 = (T_gas-1300)/1000;
        CpdT_CO_gas_2 = 35.15070*t_CO_gas_2 + (1.300095/2)*(t_CO_gas_2^2) + (-0.205921/3)*(t_CO_gas_2^3) ...
        + (0.013550/4)*(t_CO_gas_2^4) - (-3.282780)/t_CO_gas_2;
        CpdT_CO_gas = CpdT_CO_gas_1 + CpdT_CO_gas_2;
    end
    
    if T_lSc < 1300
        t_CO_lSc = (T_gas-298)/1000;
        CpdT_CO_lSc = 25.56759*t_CO_lSc + (6.096130/2)*(t_CO_lSc^2) + (4.054656/3)*(t_CO_lSc^3) ...
        + (-2.671301/4)*(t_CO_lSc^4) - (0.131021)/t_CO_lSc;
    else
        t_CO_lSc = (1300-298)/1000;
        CpdT_CO_lSc_1 = 25.56759*t_CO_lSc + (6.096130/2)*(t_CO_lSc^2) + (4.054656/3)*(t_CO_lSc^3) ...
        + (-2.671301/4)*(t_CO_lSc^4) - (0.131021)/t_CO_lSc;
        t_CO_lSc_2 = (T_lSc-1300)/1000;
        CpdT_CO_lSc_2 = 35.15070*t_CO_lSc_2 + (1.300095/2)*(t_CO_lSc_2^2) + (-0.205921/3)*(t_CO_lSc_2^3) ...
        + (0.013550/4)*(t_CO_lSc_2^4) - (-3.282780)/t_CO_lSc_2;
        CpdT_CO_lSc = CpdT_CO_lSc_1 + CpdT_CO_lSc_2;    
    end
    
    CpdT_MnO = 0.04869 * (T_lSc - 298);
    
    if T_lSc < 980
        t_Mn = (T_lSc - 298)/1000;
        CpdT_Mn = 27.24190*t_Mn + (5.237640/2)*(t_Mn^2) + (7.783160/3)*(t_Mn^3) ...
        + (-2.118501/4)*(t_Mn^4) - (-0.282113)/t_Mn;
    elseif (980 <= T_lSc) && (T_lSc < 1361)
        t_Mn = (980 - 298)/1000;
        CpdT_Mn_1 = 27.24190*t_Mn + (5.237640/2)*(t_Mn^2) + (7.783160/3)*(t_Mn^3) ...
        + (-2.118501/4)*(t_Mn^4) - (-0.282113)/t_Mn;
        t_Mn2 = (T_lSc - 980)/1000;
        CpdT_Mn_2 = 52.29870*t_Mn2 + (-28.67560/2)*(t_Mn2^2) + (21.48670/3)*(t_Mn2^3) ...
        + (-4.979850/4)*(t_Mn2^4) - (-2.432060)/t_Mn2;
        CpdT_Mn = CpdT_Mn_1 + CpdT_Mn_2;
    elseif (1361 <= T_lSc) && (T_lSc < 1412)
        t_Mn = (980 - 298)/1000;
        CpdT_Mn_1 = 27.24190*t_Mn + (5.237640/2)*(t_Mn^2) + (7.783160/3)*(t_Mn^3) ...
        + (-2.118501/4)*(t_Mn^4) - (-0.282113)/t_Mn;
        t_Mn2 = (1361 - 980)/1000;
        CpdT_Mn_2 = 52.29870*t_Mn2 + (-28.67560/2)*(t_Mn2^2) + (21.48670/3)*(t_Mn2^3) ...
        + (-4.979850/4)*(t_Mn2^4) - (-2.432060)/t_Mn2;
        t_Mn3 = (1412 - 1361)/1000;
%         CpdT_Mn_3 = 19.06450*t_Mn3 + (31.41340/2)*(t_Mn3^2) + (-14.99350/3)*(t_Mn3^3) ...
%         + (3.214741/4)*(t_Mn3^4) - (1.867090)/t_Mn3;
        CpdT_Mn_3 = 43.3 * t_Mn3;
        CpdT_Mn = CpdT_Mn_1 + CpdT_Mn_2 + CpdT_Mn_3;
    elseif (1412 <= T_lSc) && (T_lSc < 1519)
        t_Mn = (980 - 298)/1000;
        CpdT_Mn_1 = 27.24190*t_Mn + (5.237640/2)*(t_Mn^2) + (7.783160/3)*(t_Mn^3) ...
        + (-2.118501/4)*(t_Mn^4) - (-0.282113)/t_Mn;
        t_Mn2 = (1361 - 980)/1000;
        CpdT_Mn_2 = 52.29870*t_Mn2 + (-28.67560/2)*(t_Mn2^2) + (21.48670/3)*(t_Mn2^3) ...
        + (-4.979850/4)*(t_Mn2^4) - (-2.432060)/t_Mn2;
        t_Mn3 = (1412 - 1361)/1000;
%         CpdT_Mn_3 = 19.06450*t_Mn3 + (31.41340/2)*(t_Mn3^2) + (-14.99350/3)*(t_Mn3^3) ...
%         + (3.214741/4)*(t_Mn3^4) - (1.867090)/t_Mn3;
        CpdT_Mn_3 = 43.3 * t_Mn3;
        t_Mn4 = (1519 - 1412)/1000;
%         CpdT_Mn_4 = -534.1720*t_Mn4 + (679.0530/2)*(t_Mn4^2) + (-296.3700/3)*(t_Mn4^3) ...
%         + (46.42660/4)*(t_Mn4^4) - (161.3420)/t_Mn4;
        CpdT_Mn_4 = 45.8 * t_Mn4;
        CpdT_Mn = CpdT_Mn_1 + CpdT_Mn_2 + CpdT_Mn_3 + CpdT_Mn_4;
    else
        t_Mn = (980 - 298)/1000;
        CpdT_Mn_1 = 27.24190*t_Mn + (5.237640/2)*(t_Mn^2) + (7.783160/3)*(t_Mn^3) ...
        + (-2.118501/4)*(t_Mn^4) - (-0.282113)/t_Mn;
        t_Mn2 = (1361 - 980)/1000;
        CpdT_Mn_2 = 52.29870*t_Mn2 + (-28.67560/2)*(t_Mn2^2) + (21.48670/3)*(t_Mn2^3) ...
        + (-4.979850/4)*(t_Mn2^4) - (-2.432060)/t_Mn2;
        t_Mn3 = (1412 - 1361)/1000;
%         CpdT_Mn_3 = 19.06450*t_Mn3 + (31.41340/2)*(t_Mn3^2) + (-14.99350/3)*(t_Mn3^3) ...
%         + (3.214741/4)*(t_Mn3^4) - (1.867090)/t_Mn3;
        CpdT_Mn_3 = 43.3 * t_Mn3;
        t_Mn4 = (1519 - 1412)/1000;
%         CpdT_Mn_4 = -534.1720*t_Mn4 + (679.0530/2)*(t_Mn4^2) + (-296.3700/3)*(t_Mn4^3) ...
%         + (46.42660/4)*(t_Mn4^4) - (161.3420)/t_Mn4;
        CpdT_Mn_4 = 45.8 * t_Mn4;
        t_Mn5 = (T_lSc - 1519)/1000;
        CpdT_Mn_5 = 46.02400*t_Mn5 + ((1.953485e-7)/2)*(t_Mn5^2) + ((-7.567225e-8)/3)*(t_Mn5^3) ...
        + ((1.005938e-8)/4)*(t_Mn4^5) - (5.623757e-8)/t_Mn5; 
        CpdT_Mn = CpdT_Mn_1 + CpdT_Mn_2 + CpdT_Mn_3 + CpdT_Mn_4 + CpdT_Mn_5;
    end
    
    if T_lSc < 847
        t_SiO2 = (T_lSc - 298)/1000;
        CpdT_SiO2 = -6.076591*t_SiO2 + (251.6755/2)*(t_SiO2^2) + (-324.7964/3)*(t_SiO2^3) ...
        + (168.5604/4)*(t_SiO2^4) - (0.002548)/t_SiO2;
    else
        t_SiO2 = (847 - 298)/1000;
        CpdT_SiO2_1 = -6.076591*t_SiO2 + (251.6755/2)*(t_SiO2^2) + (-324.7964/3)*(t_SiO2^3) ...
        + (168.5604/4)*(t_SiO2^4) - (0.002548)/t_SiO2;
        t_SiO2_2 = (T_lSc - 847)/1000;
        CpdT_SiO2_2 = 58.75340*t_SiO2_2 + (10.27925/2)*(t_SiO2_2^2) + (-0.131384/3)*(t_SiO2_2^3) ...
        + (0.025210/4)*(t_SiO2_2^4) - (0.025601)/t_SiO2_2;
        CpdT_SiO2 = CpdT_SiO2_1 + CpdT_SiO2_2;
    end
    
    if T_lSc < 1685
        t_Si = (T_lSc - 298)/1000;
        CpdT_Si = 22.81719*(t_Si) + (3.899510/2)*(t_Si^2) + (-0.082885/3)*(t_Si^3) ...
        + (0.042111/4)*(t_Si^4) - (-0.354063)/t_Si;
    else
        t_Si = (1685 - 298)/1000;
        CpdT_Si_1 = 22.81719*(t_Si) + (3.899510/2)*(t_Si^2) + (-0.082885/3)*(t_Si^3) ...
        + (0.042111/4)*(t_Si^4) - (-0.354063)/t_Si;
        t_Si_2 = (T_lSc - 1685)/1000;
        CpdT_Si_2 = 27.19604*(t_Si_2);
        CpdT_Si = CpdT_Si_1 + CpdT_Si_2;
    end
    
    t_Cr2O3 = (T_lSc - 298)/1000;
    CpdT_Cr2O3 = 124.6550*t_Cr2O3 + (-0.337045/2)*(t_Cr2O3^2) + (5.705010/3)*(t_Cr2O3^3) ...
    + (-1.053470/4)*(t_Cr2O3^4) - (-2.030501)/t_Cr2O3;
    
    if T_lSc < 600
        t_Cr = (T_lSc - 298)/1000;
        CpdT_Cr = 7.489737*t_Cr + (71.50498/2)*(t_Cr^2) + (-91.67562/3)*(t_Cr^3) ...
        + (46.04450/4)*(t_Cr^4) - (0.138157)/t_Cr;
    else
        t_Cr = (600 - 298)/1000;
        CpdT_Cr_1 = 7.489737*t_Cr + (71.50498/2)*(t_Cr^2) + (-91.67562/3)*(t_Cr^3) ...
        + (46.04450/4)*(t_Cr^4) - (0.138157)/t_Cr;
        t_Cr2 = (T_lSc - 600)/1000;
        CpdT_Cr_2 = 18.46508*t_Cr2 + (5.477986/2)*(t_Cr2^2) + (7.904329/3)*(t_Cr2^3) ...
        + (-1.147848/4)*(t_Cr2^4) - (1.265791)/t_Cr2;
        CpdT_Cr = CpdT_Cr_1 + CpdT_Cr_2;
    end
    
    CpdT_P2O5 = 0.143 * (T_lSc - 298);
    
    if T_lSc < 1180
        CpdT_P = 0.02633 * (T_lSc -298);
    else
        CpdT_P_1 = 0.02633 * (1180 - 298);
        t_P = (T_lSc - 1180)/1000;
        CpdT_P_2 = 20.44403*t_P + (1.051745/2)*(t_P^2) + (-1.098514/3)*(t_P^3) + (0.377924/4)*(t_P^4) ...
            - (0.010645)/t_P;
        CpdT_P = CpdT_P_1 + CpdT_P_2;
    end
    
    if T_gas < 1200
        t_CO2_gas = (T_gas - 298)/1000;
        CpdT_CO2_gas = 24.99735*t_CO2_gas + (55.18696/2)*(t_CO2_gas^2) + (-33.69137/3)*(t_CO2_gas^3) ...
        + (7.948387/4)*(t_CO2_gas^4) - (-0.136638)/t_CO2_gas;
    else
        t_CO2_gas = (1200 - 298)/1000;
        CpdT_CO2_gas_1 = 24.99735*t_CO2_gas + (55.18696/2)*(t_CO2_gas^2) + (-33.69137/3)*(t_CO2_gas^3) ...
        + (7.948387/4)*(t_CO2_gas^4) - (-0.136638)/t_CO2_gas;
        t_CO2_gas_2 = (T_gas - 1200)/1000;
        CpdT_CO2_gas_2 = 58.16639*t_CO2_gas_2 + (2.720074/2)*(t_CO2_gas_2^2) + (-0.492289/3)*(t_CO2_gas_2^3) ...
        + (0.038844/4)*(t_CO2_gas_2^4) - (-6.447293)/t_CO2_gas_2;
        CpdT_CO2_gas = CpdT_CO2_gas_1 + CpdT_CO2_gas_2;
    end
    
    if T_lSc < 1200
        t_CO2_lSc = (T_lSc - 298)/1000;
        CpdT_CO2_lSc = 24.99735*t_CO2_lSc + (55.18696/2)*(t_CO2_lSc^2) + (-33.69137/3)*(t_CO2_lSc^3) ...
        + (7.948387/4)*(t_CO2_lSc^4) - (-0.136638)/t_CO2_lSc;
    else
        t_CO2_lSc = (1200 - 298)/1000;
        CpdT_CO2_lSc_1 = 24.99735*t_CO2_lSc + (55.18696/2)*(t_CO2_lSc^2) + (-33.69137/3)*(t_CO2_lSc^3) ...
        + (7.948387/4)*(t_CO2_lSc^4) - (-0.136638)/t_CO2_lSc;
        t_CO2_lSc_2 = (T_lSc - 1200)/1000;
        CpdT_CO2_lSc_2 = 58.16639*t_CO2_lSc_2 + (2.720074/2)*(t_CO2_lSc_2^2) + (-0.492289/3)*(t_CO2_lSc_2^3) ...
        + (0.038844/4)*(t_CO2_lSc_2^4) - (-6.447293)/t_CO2_lSc_2;
        CpdT_CO2_lSc = CpdT_CO2_lSc_1 + CpdT_CO2_lSc_2;
    end
    
    Cp_CH4 = 0.0586;
    CpdT_C9H20 = 0.40334 * (T_gas - 298);
    
    % ----------------- Heat of reaction -----------------
    
    % a) Fe + 1/2O2 -> FeO
    dH_Ta = r_Fe_hO2 * (dH_FeO + dH_FeS + CpdT_FeO - CpdT_Fe - 0.5*CpdT_O2_lSc);

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
    dH_Tg = r_C_hO2 * ((dH_CO-dH_CS) + CpdT_CO_lSc - CpdT_C - 0.5*CpdT_O2_lSc);

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
        CpdT_SiO2 - CpdT_Si - CpdT_O2_lSc);

    % m) 2Cr + 3/2O2 -> Cr2O3
    dH_Tm = r_2Cr_3hO2 * ((dH_Cr2O3 - 2*dH_CrS) + ...
        CpdT_Cr2O3 - 2*CpdT_Cr - 1.5*CpdT_O2_lSc);

    % Original n) removed; now n) = paper's p)
    
    % n) C9H20 + 14O2 -> 9CO2 + 10H2O
    dH_Tn = r_comb * ((9*dH_CO2 + 10*dH_H2O - dH_C9H20) + ...
        9*CpdT_CO2_gas + 10*CpdT_H2O - CpdT_C9H20 - 14*CpdT_O2_gas);

    % o) Graphite to CO2
    dH_To = (dm_el/M_C) * (dH_CO2 + CpdT_CO2_gas - CpdT_C_gas - CpdT_O2_gas);


    
    % ======================== Heat Transfer (kW) =========================
    
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
    
    VF_511 = (b1/(8*R*H1)) + (1/(2*pi))*(acos(a1/b1) - (1/(2*H1))*sqrt((((a1+2)^2)/(R^2))-4)*acos((a1*R)/b1) ...
        - (a1/(2*R*H1))*asin(R));
    VF_512 = (b2/(8*R*H2)) + (1/(2*pi))*(acos(a2/b2) - (1/(2*H2))*sqrt((((a2+2)^2)/(R^2))-4)*acos((a2*R)/b2) ...
        - (a2/(2*R*H2))*asin(R));
    
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
    FX = (bX/(8*R*X)) + (1/(2*pi))*(acos(aX/bX) - (1/(2*X))*sqrt((((aX+2)^2)/(R^2))-4)*acos((aX*R)/bX) ...
        - (aX/(2*R*X))*asin(R));
    if isnan(FX)
        FX = 0;
    end
    
    aLX = (L-X)^2 + R^2 - 1;
    bLX = (L-X)^2 - R^2 + 1;
    FLX = (bLX/(8*R*(L-X))) + (1/(2*pi))*(acos(aLX/bLX) - (1/(2*(L-X)))*sqrt((((aLX+2)^2)/(R^2))-4)*acos((aLX*R)/bLX) ...
        - (aLX/(2*R*(L-X)))*asin(R));
    
    aYXL = (Y+X-L)^2 + R^2 - 1;
    bYXL = (Y+X-L)^2 - R^2 + 1;
    FYXL = (bYXL/(8*R*(Y+X-L))) + (1/(2*pi))*(acos(aYXL/bYXL) - (1/(2*(Y+X-L)))*sqrt((((aYXL+2)^2)/(R^2))-4)*acos((aYXL*R)/bYXL) ...
        - (aYXL/(2*R*(Y+X-L)))*asin(R));
    
    aXY = (X+Y)^2 + R^2 - 1;
    bXY = (X+Y)^2 - R^2 + 1;
    FXY = (bXY/(8*R*(X+Y))) + (1/(2*pi))*(acos(aXY/bXY) - (1/(2*(X+Y)))*sqrt((((aXY+2)^2)/(R^2))-4)*acos((aXY*R)/bXY) ...
        - (aXY/(2*R*(X+Y)))*asin(R));
    
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
    J_wall = (ep2*sig*T_wall^4/1000 + (1-ep2)*(VF_21*J_roof + VF_23*J_sSc ...
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
    + dH_Ti + dH_Tj + dH_Tk + dH_Tl + dH_Tm + dH_Tn + dH_To);
    
    % Net heat flow in solid steel zone (sSc)
    % CO post combustion and Oxygen burner neglected
    Q_sSc = (Q_arc - dH_Th)*(1-K_sSclSc) + Q_lScsSc ...
    - Q_sScsSl - Q_sSclSl - Q_sScgas - Q_sScwater - Q_sScRAD;
    
    % Net heat flow in liquid metal zone (lSc)
    % CO post combustion and Oxygen burner neglected
    Q_lSc = (Q_arc - dH_Th)*K_sSclSc - Q_lScchem ...
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
    dT_lSc = Q_lSc/((m_lSc/M_Fe)*Cp_lSc);
    
    % Temperature change of sSl
    dT_sSl = (Q_sSl*(1-(T_sSl/T_melt))) / ((m_sSl/M_sSl)*Cp_sSl);
    
    % Temperature change of gas
    dT_gas = Q_gas/((m_gas/M_gas)*Cp_gas);
    
    % Temperature change of lSl
    dT_lSl = Q_lSl/((m_lSl/M_lSl)*Cp_lSl);

    % Temperature change of roof
    dT_roof = (-Q_roofRAD + (A1/(A1+A2))*Q_gaswater - phi1*Cp_H2O*(T_roof - T_water)) ...
        / (A1*d1*rho*Cp_roof);

    % Temperature change of wall
    dT_wall = (-Q_wallRAD + (A2/(A1+A2))*Q_gaswater - phi2*Cp_H2O*(T_wall - T_water)) ...
        / (A2*d2*rho*Cp_wall);
    
    T_sSc = T_sSc + dT_sSc * ts;
    T_sSl = T_sSl + dT_sSl * ts;
    T_lSc = T_lSc + dT_lSc * ts;
    T_lSl = T_lSl + dT_lSl * ts;
    T_gas = T_gas + dT_gas * ts;
    T_roof = T_roof + dT_roof * ts;
    T_wall = T_wall + dT_wall * ts;
    
    % ======================= Phase Change =======================
    
    % Injected Carbon Dissolve Rate
    dm_CL_melt = (m_CL * T_lSc * Cp_lSc * (T_air/T_melt)) / ...
        (lambda_C + Cp_C * (T_melt - T_air));
    
    % Melt rate of solid metal (kg/s)
    dm_sSc = ((Q_sSc*(T_sSc/T_melt)) / (lambda_sSc + Cp_sSc*(T_melt - T_sSc)));
    %T_lSc = (T_lSc * m_lSc + T_melt * dm_sSc * ts) / (m_lSc + dm_sSc * ts);
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
    m_Al2O3_sSc = m_Al2O3_sSc - (dm_sSc*MX_Al2O3_sSc*ts);
    m_Al2O3_lSl = m_Al2O3_lSl + (dm_sSc*MX_Al2O3_sSc*ts);
    m_CaO_sSc = m_CaO_sSc - (dm_sSc*MX_CaO_sSc*ts);
    m_CaO_lSl = m_CaO_lSl + (dm_sSc*MX_CaO_sSc*ts);
    m_SiO2_sSc = m_SiO2_sSc - (dm_sSc*MX_SiO2_sSc*ts);
    m_SiO2_lSl = m_SiO2_lSl + (dm_sSc*MX_SiO2_sSc*ts);
    m_MgO_sSc = m_MgO_sSc - (dm_sSc*MX_MgO_sSc*ts);
    m_MgO_lSl = m_MgO_lSl + (dm_sSc*MX_MgO_sSc*ts);
    m_MnO_sSc = m_MnO_sSc - (dm_sSc*MX_MnO_sSc*ts);
    m_MnO_lSl = m_MnO_lSl + (dm_sSc*MX_MnO_sSc*ts);
    m_P2O5_sSc = m_P2O5_sSc - (dm_sSc*MX_P2O5_sSc*ts);
    m_P2O5_lSl = m_P2O5_lSl + (dm_sSc*MX_P2O5_sSc*ts);
    
    % Melt rate of solid slag (kg/s)
    % Melt temperature of 1400C according to https://core.ac.uk/download/pdf/82678298.pdf
    dm_sSl = (Q_sSl*(T_sSl/1673)) / ((lambda_sSl + Cp_sSl*(1673 - T_sSl))/M_sSl);
    
    %T_lSl = (T_lSl * m_lSl + 1673 * dm_sSl * ts) / (m_lSl + dm_sSl * ts);
    
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
    m_CO = m_CO + (r_MnO_C * M_CO * ts);
    
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
    
    m_FeO_lSl = m_FeO_lSl - 1.5*(r_3FeO_2Cr * M_FeO * ts);
    m_Cr_lSc = m_Cr_lSc - (r_3FeO_2Cr * M_Cr * ts);
    m_Fe_lSc = m_Fe_lSc + 1.5*(r_3FeO_2Cr * M_Fe * ts);
    m_Cr2O3_lSl = m_Cr2O3_lSl + 0.5*(r_3FeO_2Cr * M_Cr2O3 * ts);
    
    % --------- Chromium Oxidation ---------
    
    % 2Cr + 3/2O2 -> Cr2O3
    
    m_Cr_lSc = m_Cr_lSc - (r_2Cr_3hO2 * M_Cr * ts);
    m_Cr2O3_lSl = m_Cr2O3_lSl + 0.5*(r_2Cr_3hO2 * M_Cr2O3 * ts);
    m_O2 = m_O2 - 1.5*(r_2Cr_3hO2 * M_O2 * ts);
    
    % -------- Phosphorus Oxidation --------
    
    % 5FeO + 2P -> 5Fe + P2O5
    
    m_FeO_lSl = m_FeO_lSl - 2.5*(r_5FeO_2P * M_FeO * ts);
    m_P_lSc = m_P_lSc - (r_5FeO_2P * M_P * ts);
    m_Fe_lSc = m_Fe_lSc + 2.5*(r_5FeO_2P * M_Fe * ts);
    m_P2O5_lSl = m_P2O5_lSl + 0.5*(r_5FeO_2P * M_P2O5 * ts);
    
    % ------------- Fe Oxidation -------------
    
    % Fe + 1/2O2 -> FeO
    
    m_Fe_lSc = m_Fe_lSc - (r_Fe_hO2 * M_Fe * ts);
    m_FeO_lSl = m_FeO_lSl + (r_Fe_hO2 * M_FeO * ts);
    m_O2 = m_O2 - 0.5*(r_Fe_hO2 * M_O2 * ts);
    
    % ------------- Combustion ---------------
    
    % C9H20 + 14O2 -> 9CO2 + 10H2O
    
    m_comb_sSc = m_comb_sSc - (r_comb * M_C9H20 * ts);
    m_O2 = m_O2 - 14*(r_comb * M_O2 * ts);
    m_CO2 = m_CO2 + 9*(r_comb * M_CO2 * ts);
    m_H2O = m_H2O + 10*(r_comb * M_H2O * ts);
    
    % ----------- Post Combustion ------------
    
    % CO + 1/2O2 -> CO2
    
    m_O2 = m_O2 - (r_post * M_O2) * ts;
    m_CO = m_CO - 2*(r_post * M_CO * ts);
    m_CO2 = m_CO2 + 2*(r_post * M_CO2 * ts);
    
    % --------- Electrode Oxidation ----------
    
    % C + O2 -> CO2             C here is from electrode so no eq needed
    
    m_O2 = m_O2 - ((dm_el*M_O2)/M_C) * ts;
    m_CO2 = m_CO2 + ((dm_el*M_CO2)/M_C) * ts;
    
    % ======================= Material Addition ======================
    
    % ------------ Solid Metal ------------
    
    % Total mass of solid metal
    m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_SiO2_sSc ...
        + m_Al2O3_sSc + m_CaO_sSc + m_MgO_sSc + m_MnO_sSc + m_Si_sSc + m_comb_sSc;
    
    % Addition of DRI
    m_Fe_sSc = m_Fe_sSc + (DRI_add*ts) * MX_Fe_DRI;
    m_C_sSc = m_C_sSc + (DRI_add*ts) * MX_C_DRI;
    m_SiO2_sSc = m_SiO2_sSc + (DRI_add*ts) * MX_SiO2_DRI;
    m_Al2O3_sSc = m_Al2O3_sSc + (DRI_add*ts) * MX_Al2O3_DRI;
    m_CaO_sSc = m_CaO_sSc + (DRI_add*ts) * MX_CaO_DRI;
    m_MgO_sSc = m_MgO_sSc + (DRI_add*ts) * MX_MgO_DRI;
    m_MnO_sSc = m_MnO_sSc + (DRI_add*ts) * MX_MnO_DRI;
    m_P2O5_sSc = m_P2O5_sSc + (DRI_add*ts) * MX_P2O5_DRI;
    
    % Addition of scrap
    m_Fe_sSc = m_Fe_sSc + (scr_add*ts) * MX_Fe_scr;
    m_C_sSc = m_C_sSc + (scr_add*ts) * MX_C_scr;
    m_Si_sSc = m_Si_sSc + (scr_add*ts) * MX_Si_scr;
    m_Cr_sSc = m_Cr_sSc + (scr_add*ts) * MX_Cr_scr;
    m_P_sSc = m_P_sSc + (scr_add*ts) * MX_P_scr;
    m_Mn_sSc = m_Mn_sSc + (scr_add*ts) * MX_Mn_scr;
    m_comb_sSc = m_comb_sSc + (scr_add*ts) * MX_comb_scr;
    
    T_sSc = (T_sSc * m_sSc + T_scr * scr_add * ts + T_DRI * DRI_add * ts) / (m_sSc + slg_add * ts + DRI_add * ts);
    
    % ------------ Solid Slag ------------
    m_CaO_sSl = m_CaO_sSl + (slg_add*ts) * MX_CaO_slg;
    m_MgO_sSl = m_MgO_sSl + (slg_add*ts) * MX_MgO_slg;
    m_SiO2_sSl = m_SiO2_sSl + (slg_add*ts) * MX_SiO2_slg;
    m_Al2O3_sSl = m_Al2O3_sSl + (slg_add*ts) * MX_Al2O3_slg;
    
    T_sSl = (T_sSl * m_sSl + T_slg * slg_add * ts) / (m_sSl + slg_add * ts);
    
    % ------ Extra Material Addition -----
    
    % Oxygen Lance
    m_O2 = m_O2 + O2_lance * ts;
    
    % Oxygen Post
    m_O2 = m_O2 + O2_post * ts;
    
    T_gas = (T_gas * m_gas + T_air * O2_lance * ts + T_air * O2_post * ts) / (m_gas + O2_lance * ts + O2_post * ts);
    
    % Carbon injection
    m_CL = m_CL + (C_inj * ts);
    
    % Ferro-Manganese Injection
    m_Mn_sSc = m_Mn_sSc + (FM_inj * MX_Mn_FM * ts);
    m_C_sSc = m_C_sSc + (FM_inj * MX_C_FM * ts);
    m_P_sSc = m_P_sSc + (FM_inj * MX_P_FM * ts);
    m_Si_sSc = m_Si_sSc + (FM_inj * MX_Si_FM * ts);
    m_Fe_sSc = m_Fe_sSc + (FM_inj * MX_Fe_FM * ts);
    
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
    
    dm_CO = r_FeO_CL * M_CO + r_FeO_CD * M_CO + r_C_hO2 * M_CO + r_MnO_C * M_MnO - 2*r_post * M_CO;
    dm_CO2 = r_C_O2 * M_CO2 + 9*r_comb * M_CO2 + 2*r_post * M_CO2 + dm_el* M_CO2 / M_C;
    dm_O2 = O2_lance -0.5*r_C_hO2 * M_O2 - r_C_O2 * M_O2 - r_Si_O2 * M_O2 ...
        -1.5*r_2Cr_3hO2 * M_O2 - 14*r_comb * M_O2 + O2_post - r_post * M_O2;
    dm_H2O = 10 * r_comb * M_H2O;
    
    rp = (R*T_gas/V_gas)*(dm_CO/M_CO + dm_CO2/M_CO2 + dm_O2/M_O2 + dm_H2O/M_H2O) ...
        + (R*dT_gas/V_gas)*(m_CO/M_CO + m_CO2/M_CO2 + m_O2/M_O2 + m_H2O/M_H2O);

    if mod(step, out/ts) == 0
        
        % Take out liquid metal
        lSc_out = m_lSc * 0.5;
        m_Fe_lSc = m_Fe_lSc - lSc_out * MX_Fe_lSc;
        m_C_lSc = m_C_lSc - lSc_out * MX_C_lSc;
        m_Si_lSc = m_Si_lSc - lSc_out * MX_Si_lSc;
        m_Cr_lSc = m_Cr_lSc - lSc_out * MX_Cr_lSc;
        m_Mn_lSc = m_Mn_lSc - lSc_out * MX_Mn_lSc;
        m_P_lSc = m_P_lSc - lSc_out * MX_P_lSc;

        % Take out liquid slag
        lSl_out = m_lSl * 0.5;
        m_Al2O3_lSl = m_Al2O3_lSl - lSl_out * MX_Al2O3_lSl;
        m_CaO_lSl = m_CaO_lSl - lSl_out * MX_CaO_lSl;
        m_Cr2O3_lSl = m_Cr2O3_lSl - lSl_out * MX_Cr2O3_lSl;
        m_FeO_lSl = m_FeO_lSl - lSl_out * MX_FeO_lSl;
        m_MgO_lSl = m_MgO_lSl - lSl_out * MX_MgO_lSl;
        m_MnO_lSl = m_MnO_lSl - lSl_out * MX_MnO_lSl;
        m_P2O5_lSl = m_P2O5_lSl - lSl_out * MX_P2O5_lSl;
        m_SiO2_lSl = m_SiO2_lSl - lSl_out * MX_SiO2_lSl;
        
    end 
    
    % ===================== For Graph =======================
    if mod(step, 1/ts) == 0
        gas_temp(step*ts) = T_gas;
        sSc_temp(step*ts) = T_sSc;
        sSl_temp(step*ts) = T_sSl;
        lSc_temp(step*ts) = T_lSc;
        lSl_temp(step*ts) = T_lSl;
        steel_Fe(step*ts) = MX_Fe_lSc;
        m_solid(step*ts) = m_sSc;
        m_liquid(step*ts) = m_lSc;
        m_solid_slag(step*ts) = m_sSl;
        m_liquid_slag(step*ts) = m_lSl;
        rel_pres(step*ts) = rp;
    end
    
end

kpi(1) = MX_C_lSc*100;
kpi(2) = MX_Mn_lSc*100;
kpi(3) = A_eaf * h_eafup + A_bath * h_eaflow;
kpi(4) = C_inj;
kpi(5) = FM_inj;
kpi(6) = O2_lance + O2_post;
kpi(7) = P_arc/1000;
kpi(8) = slg_add;
kpi(9) = T_lSc;
kpi(10) = (m_Fe_lSc/M_Fe)/(m_FeO_lSl/M_FeO);
kpi(11) = ((((DRI_add + slg_add + scr_add + O2_lance + O2_post + C_inj + FM_inj)*out ...
    - (m_lSc + m_lSl)/2)) * (MX_CO + MX_CO2)) / out;
kpi(12) = MX_CO2/MX_CO;
kpi(13) = m_sSc;

% Graph generation
time = linspace(1, secs, secs);

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
legend('solid', 'solid slag')
hold off

figure
plot(time, m_liquid)
hold on
plot(time, m_liquid_slag)
legend('liquid', 'liquid slag')
hold off
