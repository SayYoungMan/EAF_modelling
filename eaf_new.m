% Code based on the paper of Modeling and Validation of an Electric Arc
% Furnace Part 1 and 2 by Logar et al. 2011

clc;
clear;

% ----------------------- Variables to change -----------------------
% Carbon injection rate
C_inj = 0; % kg/s
% O2 lance rate
O2_lance = 30; % kg/s
% DRI addition rate
DRI_add = 100;
% Temperature at gas phase
T_gas = 1000; % K
% Volume of gas phase
V_gas = 10; % m^3

% ------------------------ Initial Condition ------------------------
% Mass of FeO in slag
m_FeO = 0;
% Mass of carbon in bath
m_CL = 0;
% Mass of liquid Slag
m_lSl = 0;
% Mass of liquid metal
m_lSc = 100;
% Mass of Carbon
m_C = 100;
% Mass of Si
m_Si = 200;
% Mass of SiO2
m_SiO2 = 0;
% Mass of CO
m_CO = 0;
% Mass of CO2
m_CO2 = 0;
% Mass of N2
m_N2 = 0;
% Mass of O2
m_O2 = 0;
% Temperature of liquid metal
T_lSc = 1600; % K
% Temperature of input air
T_air = 298; % K
% Melting point of steel
T_melt = 1809;


% --------------------------- Constants -----------------------------
% Heat of fusion of Carbon
lambda_C = 117; % kJ/mol
% Characteristic dimension of the duct area at the slip gap
hd = 0.65;
% Approximation of the off-gas mass flow
u1 = 15; % kg/s
% Slip-gap width
u2 = 0.3; % m
% Dimensionless constant used for improving the approximation
kU = 6.44;
% Molar amount of O2 in leak air
k_air1 = 7.3; % mol/kg
% Molar amount of N2 in leak air
k_air2 = 27.4; % mol/kg
% Ratio between the reaction rate and relative pressure
k_PR = 0.6;
% Ideal Gas Constant
R_gas = 0.008314; % kJ/mol K

% ----------------------- Reaction Rates in kg/s ---------------------
kd_CL = 15;
kd_CD = 35;
kd_C1 = 60; 
kd_C2 = 55;
kd_Si1 = 144;
kd_Si2 = 250;

k_XC = 4.9e-04;
k_XSi = 8.08e-08;

% ------------------- Specific Heat Capacities kJ/mol K -----------------
Cp_lSc = 0.047;
Cp_C = 0.03735;


% ------------------- Enthalpies of formation in kJ/mol ---------------
dH_FeO = -243;
dH_CO = -117;
dH_CO2 = -396;
dH_SiO2 = -946;

% ------------------ Fraction of the lanced oxygen used ---------------
K_O2CO = 0.05;
K_O2CO2 = 0.15;
K_O2FeO = 0.75;
K_O2SiO2 = 0.05;
K_FeODRI = 0.07;
K_FeDRI = 0.93;

% ------------------------ Molecular Weights -------------------------
% Molecular weight of Fe
M_Fe = 0.055845; % kg.mol-1
% Molecular weight of C
M_C = 0.0120107; % kg.mol-1
% Molecular weight of Si
M_Si = 0.0280855; % kg.mol-1
% Molecular weight of FeO
M_FeO = 0.071844; % kg.mol-1
% Molecular weight of SiO2
M_SiO2 = 0.06008; % kg.mol-1
% Molecular weight of O2
M_O2 = 0.032;
% Molecular weight of CO
M_CO = 0.02801;

M_lSl = 0.05;


% --------------------------- Equations -------------------------------

% Number of moles of liquid metal
XM_lSc = (m_lSc/M_Fe) + (m_C/M_C) + (m_Si/M_Si);
% Number of moles in liquid slag zones
XM_lSl = (m_lSl/M_lSl) + (m_FeO/M_FeO) + (m_SiO2/M_SiO2);

% Relative pressure in the EAF gas zone
rp = (R_gas*T_gas/V_gas)*(dm_CO/M_CO+dm_CO2/M_CO2+dm_N2/M_N2+dm_O2/M_O2) ...
    + (R_gas*dT_gas/V_gas)*(m_CO/M_CO+m_CO2/M_CO2+m_N2/M_N2+m_O2/M_O2);

% Rate of Change of C in bath
x1_d1 = C_inj;
x1_d2 = -(m_FeO*kd_CL*m_CL) / (m_lSl+m_FeO+m_SiO2);
x1_d3 = -(m_CL*T_lSc*Cp_lSc*(T_air/T_melt)) / (lambda_C+Cp_C*(T_melt-T_air));
dm_CL = x1_d1 + x1_d2 + x1_d3;

% Molecular fractions of carbon
X_C = (m_C/M_C)/XM_lSc;
X_FeO = (m_FeO/M_FeO)/XM_lSl;
Xeq_C = k_XC / X_FeO;

% Rate of change of dissolved C in the bath
x2_d1 = -kd_CD * (X_C - Xeq_C);
x2_d2 = -kd_C1 * (X_C - Xeq_C) * O2_lance * K_O2CO;
x2_d3 = -x1_d3;
x2_d5 = -kd_C2 * (X_C - Xeq_C) * O2_lance * K_O2CO2;
dm_CD = x2_d1 + x2_d2 + x2_d3 + x2_d5;

% Molecular fractions of Si
X_Si = (m_Si/M_Si)/XM_lSc;
Xeq_Si = k_XSi / X_FeO^2;

% Rate of change of Si
x3_d1 = -kd_Si1 * (X_Si - Xeq_Si);
x3_d2 = -kd_Si2 * (X_Si - Xeq_Si) * O2_lance * K_O2SiO2;
dm_Si = x3_d1 + x3_d2;

% Rate of change of FeO
x7_d1 = (M_FeO/M_C) * x1_d2;
x7_d2 = 2*(M_FeO/M_Si) * x3_d1;
x7_d3 = 2*(M_FeO/M_O2) * O2_lance * K_O2FeO;
x7_d4 = DRI_add * K_FeODRI;
x7_d5 = (M_FeO/M_C) * x2_d1;
dm_FeO = x7_d1 + x7_d2 + x7_d3 + x7_d4 + x7_d5;

% Rate of change of SiO2 in slag
dm_SiO2 = -(M_SiO2/M_Si) * dm_Si;

% Rate of change of Fe in the bath
x8_d1 = -(M_Fe/M_C) * x1_d2;
x8_d2 = -2 * (M_Fe/M_Si) * x3_d1;
x8_d3 = -(M_Fe/M_FeO) * x7_d3;
x8_d4 = DRI_add * K_FeDRI;
x8_d5 = -(M_Fe/M_C) * x2_d1;
dm_Fe = x8_d1 + x8_d2 + x8_d3 + x8_d4 + x8_d5;

% Rate of change of CO in gas zone
x9_d1 = -(hd*u1*m_CO) / ((kU*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x9_d2 = -(M_CO/M_C) * (x1_d2 + x2_d1 + x2_d2);
x9_d3 = 2*M_CO*k_air1*k_PR*rp;
x9_d4 = -(k_PR*rp*m_CO) / (m_CO+m_CO2+m_N2+m_O2);
x9_d5 = -2*(M_CO/M_O2) * O2_post * K_mCO;
dm_CO = x9_d1 + x9_d2 + x9_d4 + x9_d5;

% Rate of change of CO2 in the gas zone
x10_d1 = -(hd*u1*m_CO2) / ((kU*h2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x10_d2 = 2*M_CO2*k_air1*k_PR*rp;
x10_d3 = 2*(M_CO2/M_O2)*O2_post*K_mCO;
x10_d4 = (M_CO2/M_CH4)*CH4_inj;
x10_d5 = (M_CO2/M_C) * dm_el;
x10_d6 = -9*(M_CO2/M_C9H20)*dm_comb;
x10_d7 = -(k_PR*rp*m_CO2) / (m_CO+m_CO2+m_N2+m_O2);
x10_d8 = -(M_CO2/M_C) * x2_d5;
dm_CO2 = x10_d1 + x10_d3 + x10_d4 + x10_d5 + x10_d6 + x10_d7 + x10_d8;

% Rate of change of N2 in the gas zone
x11_d1 = -(hd*u1*m_CO2) / ((kU*h2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x11_d2 = -M_N2*k_air1*k_PR*rp;
x11_d3 = -(k_PR*rp*m_N2) / (m_CO+m_CO2+m_N2+m_O2);
dm_N2 = x11_d1 + x11_d3;

% Rate of change of oxygen
x12_d1 = -(hd*u1*m_O2) / ((kU*h2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x12_d2 = -M_O2*k_air1*k_PR*rp;
x12_d3 = -(M_O2/(2*M_CO2)) * x10_d2;
x12_d4 = O2_post * (1-K_mCO);
x12_d5 = -(M_O2/M_C) * dm_el;
x12_d6 = -14*(M_O2/M_C9H20) * dm_comb;
x12_d7 = -(k_PR*rp*m_O2) / (m_CO+m_CO2+m_N2+m_O2);
x12_d8 = -(M_O2/M_FeO) * x7_d3;
dm_O2 = x12_d1 + x12_d3 + x12_d4 + x12_d5 + x12_d6 + x12_d7 + x12_d8;

% Electrode oxidation
dm_el = 3*((R_tip*I_arc^2/3600) + R_side*A_side/3600);

% Change of combustible materials present in the solid scrap
dm_comb = -kd_comb * m_comb * (T_sSc/T_melt);

% ---------------------- Energy of Chemical Reactions --------------------
dH_ta = (x8_d3/M_Fe) * ((dH_FeO-dH_Fes) + (Cp_FeO - Cp_Fe - 0.5*Cp_O2)*(T_lSc-298));
dH_tb = ((x1_d2+x2_d2)/M_C) * ((dH_CO-dH_Cs-dH_FeO) + (Cp_Fe + Cp_CO - Cp_C - Cp_FeO)*(T_lSc-298));
dH_td = (x3_d1/M_Si) * ((dH_SiO2+dH_SiO2s-dH_FeO-dH_Sis) + (2*Cp_Fe + Cp_SiO2 - 2*Cp_FeO - 2*Cp_Si)*(T_lSc-298));
dH_tg = (x2_d2/M_C) * ((dH_CO-dH_Cs) + (Cp_CO - Cp_C - 0.5*Cp_O2)*(T_lSc-298));
dH_th = (x9_d4/M_CO) * ((dH_CO2-dH_CO) + (Cp_CO2 - Cp_CO - 0.5*Cp_O2)*(T_gas-298));
dH_ti = (x2_d5/M_C) * ((dH_CO2-dH_Cs) + (Cp_CO2 - Cp_C - Cp_O2)*(T_lSc-298));
dH_tl = (x3_d2/M_Si) * ((dH_SiO2+dH_SiO2s-dH_Sis) + (Cp_SiO2 - Cp_Si - 2*Cp_O2)*(T_lSc-298));
dH_tn = -(CH4_inj/M_CH4) * ((dH_CO2 + 2*dH_H2O - dH_CH4) + ...
    (Cp_CO2 + 2*Cp_H2O - Cp_CH4 - 2*Cp_O2));
dH_to = -(dm_el/M_C) * (dH_CO2 + (Cp_CO2 - Cp_C - Cp_O2) * (T_gas - 298));
dH_tp = -(dm_comb/M_C9H20) * ((dH_CO2 + dH_H2O - dH_C9H20) + ...
    (Cp_CO2 + Cp_H2O - Cp_C9H20 - Cp_O2)*(T_gas - 298));