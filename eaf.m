clear
clc

% ========================== Control Panel ===========================

% ---------- Time Settings ----------

% Time slice
ts = 1/1000; % 10^-3 s

% Total operating time in s
secs = 5;

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

% ---------- Other Settings ----------

% Carbon Injection Rate (kg/s)
C_inj = 0.4;

% Oxygen Lance Rate (kg/s)
O2_lance = 5;

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

% Liquid metal initial mass
m_Fe_lSc = 1000;
m_C_lSc = 20;
m_Cr_lSc = 5;
m_Mn_lSc  = 6;
m_P_lSc = 8;
m_Si_lSc = 15;

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

% Gas initial mass
m_H2O = 2;
m_O2 = 2;
m_CO = 10;
m_CO2 = 20;

% Initial mass of injected carbon
m_CL = 0;

% ------------- Initial temp. (K) --------------
T_sSc = 559.32;
T_lSc = 1850;
T_sSl = 300;
T_lSl = 1850;

% Initial Pressure
p_gas = 1.2; % atm

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
Cp_lSc = 0.047;
Cp_C = 0.02092;  % at 1800K (from NIST)

% Latent heat of fusion (kJ/mol)
lambda_C = 117;

% ------------ Constant Temperatures --------------

T_air = 298;
T_melt = 1809;

% ------------------- Others ----------------------

% Distribution of lanced oxygen
K_O2CO = 0.05;
K_O2CO2 = 0.15;
K_O2Cr2O3 = 0.015;
K_O2FeO = 0.75;
K_O2SiO2 = 0.035;

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
    r_FeO_CD = (kd_CD * (X_C_lSc - Xeq_C)) / M_C; % Logar 2011
    
    % ---------- Carbon Oxidation -----------
    
    % To carbon monoxide
    % C + 1/2 O2 -> CO
    r_C_hO2 = (kd_C1 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO) / M_C; % Logar 2011
    
    % To carbon dioxide
    % C + O2 -> CO2
    r_C_O2 = (kd_C2 * (X_C_lSc - Xeq_C) * O2_lance * K_O2CO2) / M_C; % Logar 2011
    
    % -------- MnO decarburization ----------
    
    % MnO + C -> Mn + CO
    
    % Equilibrium constant
    kX_Mn1 = 6.4 * p_CO * X_MnO_lSl;
    
    % Equilibrium mole fraction
    Xeq_MnO1 = X_Mn_lSc / kX_Mn1; % Logar 2011
    
    % Rate of reaction
    r_MnO_C = (-kd_Mn1 * (X_MnO_lSl - Xeq_MnO1)) / M_MnO;
    
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
    kX_Mn2 = 10^(2.8*B3-1.16) * ((M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2)); % Logar 2011
    
    % Equilibrium MnO mole fraction
    Xeq_MnO2 = sqrt((X_Mn_lSc^2 * X_SiO2_lSl) / (X_Si_lSc * kX_Mn2)); % Logar 2011
    
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
    
    % ========================= Phase Change =========================
    
    % Injected Carbon Dissolve Rate
    % Needs to be verified
    dm_CL_melt = (m_CL * T_lSc * Cp_lSc * (T_air/T_melt)) / ...
        (lambda_C + Cp_C * (T_melt - T_air));
    
    
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
    
end