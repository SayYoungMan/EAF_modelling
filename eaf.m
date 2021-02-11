clear
clc

% ========================== Control Panel ===========================

% ---------- Time Settings ----------

% Time slice
ts = 1/1000; % 10^-3 s

% Total operating time in s
secs = 0.1;

% ---------- DRI Settings -----------

% Temperature of DRI in K
T_DRI = 559.32;

% DRI mass addition rate in kg/s
DRI_add = 85.854;

% DRI mass fraction
MX_Fe_DRI = 0.9048;
MX_C_DRI = 0.0064;
MX_SiO2_DRI = 0.0490;
MX_Al2O3_DRI = 0.0375;
MX_CaO_DRI = 0.0014;
MX_MgO_DRI = 0.0008;
MX_MnO_DRI = 0.0001;

% --------- Scrap Settings -----------

% Temperature of scrap in K
T_scr = 559.32;

% DRI mass addition rate in kg/s
scr_add = 32.778;

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
slg_add = 10;

% Slag mass fraction
MX_CaO_slg = 0.57;
MX_MgO_slg = 0.415;
MX_SiO2_slg = 0.01;
MX_Al2O3_slg = 0.05;

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
C_inj = 0.4;

% Oxygen Lance Rate (kg/s)
O2_lance = 8;

% Power of arc (kW)
P_arc = 80000;

% EAF mass capacity (kg)
m_EAF = 200000;

% Densities (kg/m3)
rho_sSc = 900; % Kurz and Fisher 2005
rho_lSc = 7000;
rho_lSl = 3500; % Self-compacting concrete: materials properties and applications Siddique (2020) Table 10.1
rho = 7000;

% ======================= Initial Parameters =========================

% ------------- Initial mass (kg) --------------

% Solid metal initial mass
m_Fe_sSc = 950;
m_C_sSc = 10;
m_Cr_sSc = 5;
m_Mn_sSc = 5;
m_P_sSc = 5;
m_SiO2_sSc = 20;
m_Al2O3_sSc = 20;
m_CaO_sSc = 5;
m_MgO_sSc = 5;
m_MnO_sSc = 5;
m_Si_sSc = 10;
m_comb_sSc = 10;

m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_SiO2_sSc + ...
    m_Al2O3_sSc + m_CaO_sSc + m_MgO_sSc + m_MnO_sSc + m_Si_sSc + m_comb_sSc;

% Liquid metal initial mass
m_Fe_lSc = 1000;
m_C_lSc = 20;
m_Cr_lSc = 5;
m_Mn_lSc  = 6;
m_P_lSc = 8;
m_Si_lSc = 15;

m_lSc = m_Fe_lSc + m_C_lSc + m_Cr_lSc + m_Mn_lSc + m_P_lSc + m_Si_lSc;

% Solid slag initial mass
m_CaO_sSl = 20;
m_MgO_sSl = 20;
m_SiO2_sSl = 20;
m_Al2O3_sSl = 10;

% Liquid slag initial mass
m_SiO2_lSl = 10;
m_Al2O3_lSl = 20;
m_CaO_lSl = 20;
m_MgO_lSl = 20;
m_MnO_lSl = 5;
m_P2O5_lSl = 3;
m_Cr2O3_lSl = 4;
m_FeO_lSl = 20;

m_lSl = m_SiO2_lSl + m_Al2O3_lSl + m_CaO_lSl + m_MgO_lSl + m_MnO_lSl + ...
    m_P2O5_lSl + m_Cr2O3_lSl + m_FeO_lSl;

% Gas initial mass
m_H2O = 2;
m_O2 = 10;
m_CO = 10;
m_CO2 = 20;

% Initial mass of injected carbon
m_CL = 0;

% ------------- Initial temp. (K) --------------

T_sSc = 1200;
T_lSc = 1850;
T_sSl = 1000;
T_lSl = 1850;
T_gas = 1200;
T_wall = 300;
T_roof = 300;

% -------------- Initial Geometry --------------

% Cross-sectional areas
A_bath = pi * r_eafin^2;
A_eaf = pi * r_eafout^2;

% height of liquid metal
h_lSc = (m_lSc / rho_lSc) / A_bath;

% height of liquid slag
h_lSl = (m_lSl / rho_lSl) / A_bath;

% height of solid metal
V_free = (A_bath*h_eaflow) - (m_lSc/rho_lSc) - (m_lSl/rho_lSl);
if V_free >= m_sSc/rho_sSc
    h_sSc1 = 0;
    h_sSc2 = (m_sSc/rho_sSc) / A_bath;
else
    h_sSc1 = ((m_sSc/rho_sSc) - V_free) / A_eaf;
    h_sSc2 = V_free / A_bath;
end

% height of wall
h_wall = h_eafup - h_sSc1;

% Areas of roof and wall
A1 = (pi * r_eafout^2) - (pi * r_hole^2); % roof
A2 = 2 * pi * r_eafout * h_wall; % wall
A4 = pi * r_eafin^2;

% Cone geometry
d_coneout = 0.1;
d_conein = d_coneout/2;
h_cone = (sqrt(2)/2) * d_coneout;

% Arc height
h_arc = h_eafup - h_electrode - (h_sSc2 - h_cone);

% ----------------- Others ---------------------

% Initial Pressure
p_gas = 1.2; % atm

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
M_gas = 0.35;
M_Al2O3 = 0.10196;
M_Al = 0.02698;
M_H2O = 0.01802;
M_sSl = 0.0606;
M_lSl = 0.0606;

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

% Cooling water flowrate
phi1 = 10;
phi2 = 5;

for step = 1:secs/ts
    
    % ====================== Mole, Mass Fraction =====================
    
    % ------------ Solid Metal ------------
    
    % Total mass of solid metal
    m_sSc = m_Fe_sSc + m_C_sSc + m_Cr_sSc + m_Mn_sSc + m_P_sSc + m_SiO2_sSc ...
        + m_Al2O3_sSc + m_CaO_sSc + m_MgO_sSc + m_MnO_sSc + m_Si_sSc + m_comb_sSc;
    
    % Total mole of solid metal
    XM_sSc = (m_Fe_sSc/M_Fe) + (m_C_sSc/M_C) + (m_Cr_sSc/M_Cr) + (m_Mn_sSc/M_Mn) ... 
        + (m_P_sSc/M_P) + (m_SiO2_sSc/M_SiO2) + (m_Al2O3_sSc/M_Al2O3) + ...
        (m_CaO_sSc/M_CaO) + (m_MgO_sSc/M_MgO) + (m_MnO_sSc/M_MnO) + ...
        (m_Si_sSc/M_Si) + (m_comb_sSc/M_C9H20);
    
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
    
    % Equilibrium constant
    kX_Mn1 = 6.4 * p_CO * X_MnO_lSl;
    
    % Equilibrium mole fraction
    Xeq_MnO1 = X_Mn_lSc / kX_Mn1; % Logar 2012
    
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
    
    % Equilibrium fraction
    MXeq_Si = a_SiO2_bas / (K_Si * MX_FeO_lSl^2); % Turkogan 1996
    Xeq_Si = ((MXeq_Si * m_lSc) / M_Si) / XM_lSc; 
    % Xeq_Si = 8.08e-08* ((m_lSl*M_FeO)/(m_FeO_lSl*M_lSl) + (m_SiO2_lSl*M_FeO) ...
    %     /(m_FeO_lSl*M_SiO2) + 1)^2; % Bekkar 1999

    % Rate of reaction
    r_2FeO_Si = (kd_Si1 * (X_Si_lSc - Xeq_Si)) / M_Si;
    
    % --------- Silicon Oxidation ---------
    
    % Si + O2 -> SiO2
    
    r_Si_O2 = (kd_Si2 * (X_Si_lSc - Xeq_Si) * O2_lance * K_O2SiO2) / M_Si;
    
    % ------- Si reaction with MnO --------
    
    % 2MnO + Si -> 2Mn + SiO2
    
    % Reaction equilibrium constant
    kX_Mn2 = 10^(2.8*B3-1.16) * ((M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2)); % Logar 2012
    
    % Equilibrium MnO mole fraction
    Xeq_MnO2 = sqrt((X_Mn_lSc^2 * X_SiO2_lSl) / (X_Si_lSc * kX_Mn2)); % Logar 2012
    
    % Rate of reaction
    r_2MnO_Si = (kd_Mn2 * (X_MnO_lSl - Xeq_MnO2)) / M_MnO;
    
    % ------- Mn reaction with FeO ---------
    
    % Mn + FeO -> MnO + Fe
    
    % Equilibrium constant
    K_FeMn = 1.8; % Turkdogan 1996 (Given that B = 2.5-4.0 and 1600 - 1650 C)
    kX_Mn = K_FeMn * (M_FeO * M_Mn * 100) / (M_MnO * M_Fe);
    
    % Equilibrium fraction
    Xeq_Mn = X_MnO_lSl / (X_FeO_lSl * kX_Mn);
    
    % Rate of reaction
    r_FeO_Mn = (kd_Mn * (X_Mn_lSc - Xeq_Mn)) / M_Mn;
    
    % ------- Cr reaction with FeO ---------
    
    % 3FeO + 2Cr -> 3Fe + Cr2O3
    
    % Equilibrium constants
    K_FeCr = 0.3; % Trukdogan
    kX_Cr = K_FeCr * (M_Cr*M_FeO*100) / (M_Cr2O3*M_Fe);
    
    % Equilibrium mole fraction
    Xeq_Cr = X_Cr2O3_lSl / (X_FeO_lSl*kX_Cr);
    
    % Rate of reaction
    r_3FeO_2Cr = (kd_Cr1 * (X_Cr_lSc - Xeq_Cr)) / M_Cr;
    
    % --------- Chromium Oxidation ---------
    
    % 2Cr + 3/2O2 -> Cr2O3
    
    r_2Cr_3hO2 = (2*kd_Cr2 * (X_Cr_lSc - Xeq_Cr) * O2_lance * K_O2Cr2O3) / M_Cr;
    
    % ----- Phosphorus reaction with FeO -----
    
    % 5FeO + 2P -> 5Fe + P2O5
    
    partition = 10 ^ (1.97*X_CaO_lSl + 2.0*X_FeO_lSl - 2.04*X_SiO2_lSl + 6713/T_lSl - 1.84); % Basu, 2007
    meq_P = (2 * (m_P2O5_lSl/M_P2O5) * M_P) / partition;
    Xeq_P = (meq_P/M_P) / XM_lSc;
    r_5FeO_2P = (kd_P * (X_P_lSc - Xeq_P)) / M_P;
    
    % ------------- Fe Oxidation -------------
    
    % Fe + 1/2O2 -> FeO
    
    r_Fe_hO2 = (2 * O2_lance * K_O2FeO) / M_O2;
    
    % ------------- Combustion ---------------
    
    % C9H20 + 14O2 -> 9CO2 + 10H2O
    
    r_comb = (kd_comb * m_comb_sSc * (T_sSc/T_melt)) / M_C9H20;
    
    % =================== Reaction Heat Transfer ====================
    
    % ------------ Dynamic heat capacity calculation ----------
    % Units in kJ/mol K
    
    Cp_H2O = 0.050;
    Cp_C = 0.02093;
    Cp_FeO = 0.0682;
    Cp_Fe = 0.04257;
    Cp_O2 = 0.0370;
    Cp_CO = 0.0350;
    Cp_MnO = 0.04869;
    Cp_Mn = 0.04602;
    Cp_SiO2 = 0.06877;
    Cp_Si = 0.0248;
    Cp_Cr2O3 = 0.0831;
    Cp_Cr = 0.02836;
    Cp_P2O5 = 0.143;
    Cp_P = 0.02633;
    Cp_CO2 = 0.04381;
    Cp_CH4 = 0.0586;
    Cp_C9H20 = 0.40334;
    
    % a) Fe + 1/2O2 -> FeO
    dH_Ta = r_Fe_hO2 * (dH_FeO + (Cp_FeO - Cp_Fe - 0.5*Cp_O2) * (T_lSc - 298));

    % b) FeO + C -> Fe + CO
    dH_Tb = (r_FeO_CL + r_FeO_CD) * (dH_CO - dH_CS - dH_FeO + (Cp_Fe + Cp_CO - Cp_C - Cp_FeO) * (T_lSc - 298));

    % c) FeO + Mn -> Fe + MnO
    dH_Tc = r_FeO_Mn * (dH_MnO - dH_FeO - dH_MnS + (Cp_Fe + Cp_MnO - Cp_FeO - Cp_Mn) * (T_lSc - 298));

    % d) 2FeO + Si -> 2Fe + SiO2
    dH_Td = r_2FeO_Si * ((dH_SiO2+dH_SiO2S-2*dH_FeO-dH_SiS) + ...
        (2*Cp_Fe + Cp_SiO2 - 2*Cp_FeO - Cp_Si)*(T_lSc-298));

    % e) 3FeO + 2Cr -> 3Fe + Cr2O3
    dH_Te = r_3FeO_2Cr * (dH_Cr2O3 - 3*dH_FeO - 2*dH_CrS + ...
        (3*Cp_Fe + Cp_Cr2O3 - 3*Cp_FeO - 2*Cp_Cr) * (T_lSc-298));

    % f) 5FeO + 2P -> 5Fe + P2O5
    dH_Tf = r_5FeO_2P * (dH_P2O5 - 5*dH_FeO - 2*dH_PS + ...
        (5*Cp_Fe + Cp_P2O5 - 5*Cp_FeO - 2*Cp_P) * (T_lSc-298));

    % g) C + 1/2O2 -> CO
    dH_Tg = r_C_hO2 * ((dH_CO-dH_CS) + (Cp_CO - Cp_C - 0.5*Cp_O2)*(T_lSc-298));

%     % h) CO + 1/2O2 -> CO2
%     dH_Th = r_CO_hO2 * ((dH_CO2-dH_CO) + (Cp_CO2 - Cp_CO - 0.5*Cp_O2)*(T_gas-298));

    % i) C + O2 -> CO2
    dH_Ti = r_C_O2 * ((dH_CO2-dH_CS) + (Cp_CO2 - Cp_C - Cp_O2)*(T_lSc-298));

    % j) MnO + C -> Mn + CO
    dH_Tj = r_MnO_C * ((dH_CO + dH_MnS - dH_MnO - dH_CS) + ...
        (Cp_Mn + Cp_CO - Cp_MnO - Cp_C) * (T_lSc - 298));

    % k) 2MnO + Si -> 2Mn + SiO2
    dH_Tk = r_2MnO_Si * ((dH_SiO2 + dH_SiO2S + 2*dH_MnS - 2*dH_MnO - dH_SiS) + ...
        (2*Cp_Mn + Cp_SiO2 - Cp_Si - 2*Cp_MnO) * (T_lSc - 298));

    % l) Si + O2 -> SiO2
    dH_Tl = r_Si_O2 * ((dH_SiO2 + dH_SiO2S - dH_SiS) + ...
        (Cp_SiO2 - Cp_Si - 2*Cp_O2) * (T_lSc - 298));

    % m) 2Cr + 3/2O2 -> Cr2O3
    dH_Tm = r_2Cr_3hO2 * ((dH_Cr2O3 - 2*dH_CrS) + ...
        (Cp_Cr2O3 - 2*Cp_Cr - 1.5*Cp_O2) * (T_lSc - 298));

    % n) CH4 + 2O2 -> CO2 + 2H2O
%     dH_Tn = -(CH4_inj/M_CH4) * ((dH_CO2 + 2*dH_H2O - dH_CH4) + ...
%         (Cp_CO2 + 2*Cp_H2O - Cp_CH4 - 2*Cp_O2) * (T_gas - 298));

%     % o) Graphite to CO2
%     dH_To = -(dm_el/M_C) * (dH_CO2 + (Cp_CO2 - Cp_C - Cp_O2) * (T_gas - 298));

    % p) C9H20 + 14O2 -> 9CO2 + 10H2O
    dH_Tp = r_comb * ((9*dH_CO2 + 10*dH_H2O - dH_C9H20) + ...
        (9*Cp_CO2 + 10*Cp_H2O - Cp_C9H20 - 14*Cp_O2) * (T_gas - 298));
    
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
    
    % ====================== Reactor Geometry ====================
    
    % Areas of roof and wall
    A1 = (pi * r_eafout^2) - (pi * r_hole^2); % roof
    A2 = 2 * pi * r_eafout * h_wall; % wall
    
    % Surface area of sSc and lSc
    A3 = (pi * r_eafout^2) - (pi * (d_coneout/2)^2) + (pi*0.75*d_coneout*sqrt(h_cone + d_coneout/4));
    A4 = pi * r_eafin^2;
    
    
    
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
    a2 = H1^2 + R^2 - 1;
    b1 = H1^2 - R^2 + 1;
    b2 = H1^2 - R^2 + 1;
    
    VF_511 = (b1/(8*R*H1)) + (1/(2*pi))*(acos(a1/b1) - (1/(2*H1))*sqrt((a1+2)^2/(R^2)-4) ...
        *acos((a1*R)/b1) - (a1/(2*R*H1))*asin(R));
    VF_512 = (b1/(8*R*H2)) + (1/(2*pi))*(acos(a2/b2) - (1/(2*H2))*sqrt((a2+2)^2/(R^2)-4) ...
        *acos((a2*R)/b2) - (a2/(2*R*H2))*asin(R));
    
    A511 = 2*pi*r_electrode*h_wall;
    A512 = 2*pi*r_electrode*h_electrode;
    A513 = 2*pi*r_electrode*h_arc;
    
    VF_51 = (VF_511*A511 - VF_512*A512) / A513;
    
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
    
    VF_52 = (X/L)*FX + ((L-X)/L)*(1-FLX) + ((Y+X-L)/L)*FYXL - ((X+Y)/L)*FXY;
    
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
    H = h_wall / r_hole;
    R2 = r_eafout / r_hole;
    R3 = d_coneout/2 / r_hole;
    R4 = d_conein/2 / r_hole;
    
    % TODO: This equation needs checking especially the value of R2 and R4
    VF_131 = 1/(2*(R2^2-1)) * (sqrt((R2^2+R3^2+H^2)^2 - (2*R3*R2)^2) - ...
        sqrt((R2^2+R4^2+H^2)^2 - (2*R2*R4)^2) + sqrt((1+R4^2+H^2)^2 - (2*R4^2)^2) ...
        - sqrt((1+R3^2+H^2)^2 - (2*R3^2)^2));
    
    H = h_cone / (d_conein/2);
    R = (d_coneout/2) / (d_conein/2);
    X = 1 + R^2 + H^2;
    
    VF_132 = (2*R^2-X+sqrt(X^2-4*R^2))/(2*sqrt(X-2*R)*(1+R));
    
    VF_13 = VF_131 + VF_132 * VF_41;
    
    % VF_31 sSc -> Roof
    VF_31 = VF_13 * (A1/A3);
    
    % VF_32 sSc -> Wall
    R = r_eafout / (d_coneout/2);
    H = h_wall / (d_coneout/2);
    
    VF_321 = VF_132 * (1/2) * (1-R^2-H^2+sqrt((1+R^2+H^2)^2-4*R^2));
    VF_322 = (1/2)*(1 + (1/(R^2-1)) * (H*sqrt(4*R^2+H^2) - ...
        sqrt((1+R^2+H^2)^2 - 4*R^2)));
    
    VF_32 = VF_321 + VF_322;
    
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
    
    % Neglected
    VF_12 = 0; % small due to low temp.
    VF_21 = 0;
    VF_34 = 0; % dominated by conduction
    VF_43 = 0;
    
    % ========================= Radiosity ========================
    Q_arcRAD = 0.75 * P_arc;
    
    % Radiosity of roof
    J_roof = (ep1*sig*T_roof^4/1000 - (1-ep1)*(VF_12*J_wall + VF_13*J_sSc ...
        + VF_14*J_lSc + VF_15*Q_arcRAD));

    % Radiosity of wall
    J_wall = (ep2*sig*T_wall^4/1000 - (1-ep2)*(VF_21*J_wall + VF_23*J_sSc ...
        + VF_24*J_lSc + VF_25*Q_arcRAD));

    % Radiosity of sSc
    J_sSc = (ep3*sig*T_sSc^4/1000 - (1-ep3)*(VF_31*J_roof + VF_32*J_wall ...
        + VF_35*Q_arcRAD));

    % Radiosity of lSc
    J_lSc = (ep4*sig*T_lSc^4/1000 - (1-ep4)*(VF_41*J_roof + VF_42*J_wall ...
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
    
    Q_lScchem = dH_Ta + dH_Tb + dH_Tc + dH_Td + dH_Te + dH_Tf + dH_Tg ...
    + dH_Ti + dH_Tj + dH_Tk + dH_Tl + dH_Tm + dH_Tp;
    
    % Net heat flow in solid steel zone (sSc)
    % CO post combustion and Oxygen burner neglected
    Q_sSc = (Q_arc)*(1-K_sSclSc) + Q_lScsSc ...
        - Q_sScsSl - Q_sSclSl - Q_sScgas - Q_sScwater - Q_sScRAD;
    
    % Net heat flow in liquid metal zone (lSc)
    % CO post combustion and Oxygen burner neglected
    Q_lSc = (Q_arc)*K_sSclSc + Q_lScchem ...
        - Q_lScsSc - Q_lScsSl - Q_lSclSl - Q_lScgas - Q_lScwater - Q_lScRAD;
    
    % Net heat flow in solid slag zone
    Q_sSl = Q_sScsSl + Q_lScsSl - Q_sSlwater;
    
    % Net heat flow in liquid slag zone
    Q_lSl = Q_lSclSl + Q_sSclSl - Q_lSlgas - Q_lSlwater;
    
    % Gas zone energy balance
    Q_gas = Q_arcgas + + Q_sScgas + Q_lScgas + Q_lSlgas - Q_gaswater;
    
    
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
    
    % ======================= Phase Change =======================
    
    % Injected Carbon Dissolve Rate
    % Needs to be verified
    dm_CL_melt = (m_CL * T_lSc * Cp_lSc * (T_air/T_melt)) / ...
        (lambda_C + Cp_C * (T_melt - T_air));
    
    % Melt rate of solid metal (kg/s)
    dm_sSc = -((Q_sSc*(T_sSc/T_melt)) / (lambda_sSc + Cp_sSc*(T_melt - T_sSc))) * M_Fe;
    
    % Melt rate of solid slag (kg/s)
    dm_sSl = -(Q_sSl*(T_sSl/T_melt)) / ((lambda_sSl + Cp_sSl*(T_melt - T_sSl))/M_sSl);
    
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
    
    % To carbon dioxide
    % C + O2 -> CO2
    m_C_lSc = m_C_lSc - (r_C_O2 * M_C * ts);
    m_CO2 = m_CO2 + (r_C_O2 * M_CO2 * ts);
    
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
    m_CO = m_CO + 9*(r_comb * M_CO2 * ts);
    m_H2O = m_H2O + 10*(r_comb * M_H2O * ts);
    
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
    
    % Addition of scrap
    m_Fe_sSc = m_Fe_sSc + (scr_add*ts) * MX_Fe_scr;
    m_C_sSc = m_C_sSc + (scr_add*ts) * MX_C_scr;
    m_Si_sSc = m_Si_sSc + (scr_add*ts) * MX_Si_scr;
    m_Cr_sSc = m_Cr_sSc + (scr_add*ts) * MX_Cr_scr;
    m_P_sSc = m_P_sSc + (scr_add*ts) * MX_P_scr;
    m_Mn_sSc = m_Mn_sSc + (scr_add*ts) * MX_Mn_scr;
    m_comb_sSc = m_comb_sSc + (scr_add*ts) * MX_comb_scr;
    
    % ------------ Solid Slag ------------
    m_CaO_sSl = m_CaO_sSl + (slg_add*ts) * MX_CaO_slg;
    m_MgO_sSl = m_MgO_sSl + (slg_add*ts) * MX_MgO_slg;
    m_SiO2_sSl = m_SiO2_sSl + (slg_add*ts) * MX_SiO2_slg;
    m_Al2O3_sSl = m_Al2O3_sSl + (slg_add*ts) * MX_Al2O3_slg;
    
    % ------ Extra Material Addition -----
    
    % Carbon injection
    m_CL = m_CL + (C_inj * ts);
    
    % ====================== For Dev =======================
    
    O2_lance_use = ((r_C_hO2/2) + r_C_O2 + 1.5*r_2Cr_3hO2 + (r_Fe_hO2/2) ...
        + r_Si_O2 + 14*r_comb) * M_O2;
    
end