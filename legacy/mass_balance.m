clc
clear

sec = 10;
tc = 1/1000;

C_inj = 0.1;
V_gas = 45;

% Initial Mass (kg)
m_sSc = 32500;
m_lSc = 7500;
m_scrap = m_sSc + m_lSc;

m_lSl = 800;
m_sSl = 200;
m_slag = m_lSl + m_sSl;

m_gas = V_gas * 1.225;

m_C = m_scrap * 0.004;
m_CL = 0;
m_Fe = m_scrap * 0.9705;
m_Si = m_scrap * 0.006;
m_Cr = m_scrap * 0.002;
m_Mn = m_scrap * 0.006;
m_P = m_scrap * 0.0005;
m_FeO = m_Fe * 0.001;
m_SiO2 = m_slag * 0.007;
m_MnO = 1;
m_Cr2O3 = 1;
m_P2O5 = 1;
m_Al2O3 = m_slag * 0.0045;

m_CaO = m_slag * 0.567;
m_MgO = m_slag * 0.412;

m_CO = m_gas*0.005;
m_CO2 = m_gas*0.005;
m_N2 = m_gas*0.78;
m_O2 = m_gas*0.21;

m_comb = m_scrap * 0.011;

% Reaction rates in kg/s
kd_CL = 15;
kd_CD = 35;
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

% Initial Temperature
T_sSc = 1050;
T_lSc = 1800;
T_sSl = 950;
T_lSl = 1800;
T_gas = 2000;
T_wall = 308;
T_roof = 321;

% Important temperatures K
T_melt = 1809;
T_water = 298;
T_air = 298;

% ------------------------ Molecular Weights -------------------------
% kg.mol-1
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
M_MnO = 0.079374;
M_Cr2O3 = 0.15199;
M_P2O5 = 0.283886;
M_CaO = 0.05608;
M_MgO = 0.04304;
M_C9H20 = 0.1282;
M_N2 = 0.028013;
M_gas = 0.35;
M_Al2O3 = 0.10196;

M_sSl = 0.056;
M_lSl = 0.056;

% Specific heat capacity kJ/molK
Cp_sSc = 0.039;
Cp_lSc = 0.047;
Cp_sSl = 0.025;
Cp_lSl = 0.047;
Cp_gas = 0.030;
Cp_H2O = 0.074;
Cp_C = 0.03735;
Cp_FeO = 0.060;
Cp_Fe = 0.0444;
Cp_O2 = 0.0344;
Cp_CO = 0.0305;
Cp_MnO = 0.04869;
Cp_Mn = 0.02637;
Cp_SiO2 = 0.06877;
Cp_Si = 0.0248;
Cp_Cr2O3 = 0.0831;
Cp_Cr = 0.02836;
Cp_P2O5 = 0.143;
Cp_P = 0.02633;
Cp_CO2 = 0.04381;
Cp_CH4 = 0.0586;
Cp_C9H20 = 0.40334;

Cp_roof = 0.65; % kJ/kg K
Cp_wall = 0.96; % kJ/kg K

% Latent heat of fusion
lambda_sSc = 15.4; % kJ/mol
lambda_sSl = 12.66; % kJ/mol
lambda_C = 117; % kJ/mol

for step = 1:sec/tc
    
% Number of moles of liquid metal (mol)
XM_lSc = (m_lSc/M_Fe) + (m_C/M_C) + (m_Si/M_Si) + (m_Cr/M_Cr) + ...
    (m_Mn/M_Mn) + (m_P / M_P);
% Number of moles in liquid slag zones (mol)
XM_lSl = (m_lSl/M_lSl) + (m_FeO/M_FeO) + (m_SiO2/M_SiO2) + (m_MnO/M_MnO) ...
    + (m_Cr2O3/M_Cr2O3) + (m_P2O5/M_P2O5) + (m_MgO/M_MgO) + (m_CaO/M_CaO) + (m_Al2O3/M_Al2O3);
% Number of moles in gas zone (mol)
XM_gas = (m_O2/M_O2) + (m_CO/M_CO) + (m_CO2/M_CO2) + (m_N2/M_N2);

% Mole fractions
X_C = (m_C/M_C) / XM_lSc;
X_Si = (m_Si/M_Si) / XM_lSc;
X_lSc = (m_lSc/M_Fe) / XM_lSc;
X_FeO = (m_FeO/M_FeO) / XM_lSl;
X_CaO = (m_CaO/M_CaO) / XM_lSl;
X_SiO2 = (m_SiO2/M_SiO2) / XM_lSl;
X_MgO = (m_MgO/M_MgO) / XM_lSl;
X_Mn = (m_Mn/M_Mn) / XM_lSc;
X_MnO = (m_MnO/M_MnO) / XM_lSl;
X_Cr2O3 = (m_Cr2O3/M_Cr2O3) / XM_lSl;
X_P2O5 = (m_P2O5/M_P2O5) / XM_lSl;
X_Cr = (m_Cr/M_Cr) / XM_lSc;
X_P = (m_P/M_P) / XM_lSc;
X_CO = (m_CO/M_CO) / XM_gas;
X_N2 = (m_N2/M_N2) / XM_gas;
X_O2 = (m_O2/M_O2) / XM_gas;
X_CO2 = (m_CO2/M_CO2) / XM_gas;

% Equilibrium Molar Fraction

    
% Rate of change of carbon injected
x1_d1 = C_inj;
x1_d2 = -(m_FeO * kd_CL * m_CL) / ...
    (m_lSl + m_FeO + m_SiO2 + m_MnO + M_Cr2O3 + m_P2O5 + m_Al2O3);
x1_d3 = -(m_CL * T_lSc * Cp_lSc * (T_air/T_melt)) / (lambda_C + Cp_C * (T_melt - T_air));
dm_CL = x1_d1 + x1_d2 + x1_d3;
m_CL = m_CL + dm_CL * tc;

% Rate of change of carbon dissolved
x2_d1 = -kd_CD * (X_C - Xeq_C);
end