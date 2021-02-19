clear
clc

% ---------------------- Initial Conditions -----------------------
% Time slice
tc = 1/1000; % 10^-3 s
secs = 0.001; % Total operating time in s

% Arc Current
I_arc = 44; % kA

% Gas zone volume
V_gas = 45; % m^3

% Initial Mass (kg)
m_scrap = 25000;
m_DRI = 95000;
solid_ratio = 0.1;
m_slag = 5000;

% Initial Mass (kg)
m_sSc = (m_scrap + m_DRI);
m_lSc = 1;
m_lSl = m_slag * 0.2;
m_sSl = m_slag * 0.8;
m_gas = V_gas * 1.225;
m_C = m_scrap * 0.004 + m_DRI * 0.0222;
m_CL = 0;
m_Fe = m_scrap * 0.9705 + m_DRI * 0.9;
m_Si = m_scrap * 0.006;
m_Cr = m_scrap * 0.002;
m_Mn = m_scrap * 0.006;
m_P = m_scrap * 0.0005;
m_FeO = m_Fe * 0.001;
m_SiO2 = m_slag * 0.007 + m_DRI * 0.05386;
m_MnO = 0.1;
m_Cr2O3 = 0.1;
m_P2O5 = 1;
m_CaO = m_slag * 0.567;
m_MgO = m_slag * 0.412;
m_CO = m_gas*0.05;
m_CO2 = m_gas*0.05;
m_N2 = m_gas*0.78;
m_O2 = m_gas*0.21;

m_comb = m_scrap * 0.001;

% C Injection rate kg/s
C_inj = 5;

% Oxygen burner rate kg/s
CH4_inj = 0.1;

% O2 lance rate kg/s
O2_lance = 1;

% O2 rate for CO post-combustion kg/s
O2_post = 0.4;
K_mCO = m_CO / (m_CO+m_CO2);

% DRI Addition rate kg/s
DRI_add = 0;

% Initial Temperature
T_sSc = 1231.27;
T_lSc = 1800;
T_sSl = 1245;
T_lSl = 1800;
T_gas = 800;
T_wall = 300;
T_roof = 300;

% Relative Pressure
rp = 1;

% Initial radiative heat flows
Q_sScRAD = 0; %kW 
Q_lScRAD = 0; %kW

% Initial radiosity
J_roof = 0;
J_wall = 0;
J_sSc = 0;
J_lSc = 0;

% ------------------------- Variables -----------------------------
% EAF mass capacity
m_EAF = 105000; % kg
K_post = 0.7;

% Arc Power kW
P_arc = 30000;

% ------------------------- Constants -----------------------------
% Dimensionless constant for approximation
k_U = 6.44;

% Ratio between reaction rate and relative pressure
k_PR = 0.6;

% Universal Gas Constant
R = 8.314; % J/mol K

% Stefan Boltzmann Constant
sig = 5.67e-08; % W.m-2.K-4s-1

% Emissivity
ep1 = 0.85; % Roof
ep2 = 0.85; % Wall
ep3 = 0.80; % sSc
ep4 = 0.40; % lSc

% Burner efficiency
K_burn = 0.7;

% Material density of roof and wall
rho = 7000; % kg/m^3

% Thickness of roof and wall (m)
d1 = 0.3;
d2 = 0.45;

% Water flow rates kg/s
phi1 = 60;
phi2 = 130;

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

% Molar amount of O2 and N2 in leak air
k_air1 = 7.3; % mol/kg
k_air2 = 27.4; % mol/kg

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

% Enthalpies of formation kJ/mol
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

% Important temperatures K
T_melt = 1809;
T_water = 298;
T_air = 298;

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

% Fractions of lanced oxygen used
K_O2CO = 0.05;
K_O2CO2 = 0.15;
K_O2Cr2O3 = 0.015;
K_O2FeO = 0.75;
K_O2SiO2 = 0.035;

% Fraction in DRI stream
K_FeODRI = 0.07;
K_FeDRI = 0.93;

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

M_sSl = 0.056;
M_lSl = 0.056;

% ------------------------ Reactor Geometry -----------------------
% Radius
r_eafout = 6.6/2; % m
r_hole = 3.4/2; % m
r_eafin = 4.9/2;

R_tip = 0.02; % kg/kA hr
R_side = 10; % kg/m^2 hr

% Height
h_eafup = 2.9;
h_eaflow = 1;
h_wall = 0.825;
h_arc = 0.6; % Obtained from other research
h_bath = h_eaflow/2;
h_scrap = h_eafup + h_eaflow - h_wall;
h_cone = h_arc;

% Area
A1 = (r_eafout^2 - r_hole^2) * pi;
A2 = pi*r_eafout*2*h_wall;
A3 = (r_eafout^2 - r_hole^2) * pi;
A4 = pi*r_eafin^2;
A5 = pi*r_hole*2*h_arc;
A_side = 0.35;

% Characteristic dimension of the duct area at the slip gap
hd = 0;
% Approximation of off-gas mass flowrate
u1 = 15; % kg/s
% Slip-gap width
u2 = 0.3; % m

% ---------------------- View Factor Equations --------------------
K_sSclSc = 0.5 * tanh(5*(h_bath-h_scrap+h_cone)) + 0.5;
ratio = (r_eafout^2*pi) / (r_eafin^2*pi);
r1 = r_eafout/h_wall;
rho1 = (sqrt(4*r1^2 + 1) - 1) / r1;

L = h_eafup + h_eaflow;
R1 = r_eafout/L;
R4 = r_eafin/L;
S1 = 1 + (1+R4^2)/(R1^2);

VF12 = (rho1/(2*r1))/2;
VF13 = (1 - (rho1/(2*r1)))/2;
VF14 = (1/2)*(S1 - sqrt(S1^2 - 4*(r_eafout/r_eafin)^2))/2;
VF15 = 1 - VF12 - VF13 - VF14;

R2 = r_eafout/r_eafin;
H2 = h_wall/r_eafin;


VF21 = (rho1/4)/ratio;
VF22 = 1 - (rho1/2);
VF23 = (rho1/4)/ratio;
VF42 = (1/2) * (1-R2^2-H2^2 + sqrt((1+R2^2+H2^2)^2 - 4*R2^2));
VF24 = (A4/A2)*VF42;
VF25 = 1 - VF21 - VF22 - VF23 - VF24;

L3 = h_eaflow;
R3 = r_eafout/L3;
R34 = r_eafin/L3;
S3 = 1 + (1+R34^2)/(R3^2);

% VF51 = (A1/A5) * VF15;
% VF52 = (A2/A5) * VF25;
% VF53 = (1-VF51-VF52) * K_sSclSc;
% VF54 = (1-VF51-VF52) * (1-K_sSclSc);

VF31 = VF13;
VF32 = VF12;
VF34 = (1/2)*(S3 - sqrt(S3^2 - 4*(r_eafout/ratio/r_eafin)^2))/ratio;
VF35 = 1 - VF31 - VF32 - VF34;

VF41 = (A1/A4) * VF14;
% VF43 = (A3/A4) * VF34;
% VF45 = (A5/A4) * VF54;

% Custom definition for View factor
% VF12 = 0;
% VF13 = 0;
% VF14 = 0;
% VF15 = 0;
% VF21 = 0;
% VF22 = 0;
% VF23 = 0;
% VF24 = 0;
% VF25 = 0;
% VF31 = 0;
% VF32 = 0;
% VF34 = 0;
% VF35 = 0;
% VF41 = 0;
% VF42 = 0;
VF45 = 0.2;
VF51 = 0.1;
VF52 = 0.2;
VF53 = 0.3;
VF54 = 0.4;


% ------------------------- Arrays for graph ------------------------
% gas_temp = zeros(1,secs);
% sSc_temp = zeros(1,secs);
% sSl_temp = zeros(1,secs);
% lSc_temp = zeros(1,secs);
% lSl_temp = zeros(1,secs);
% steel_Fe = zeros(1,secs);

for step = 1:secs/tc
% ----------------------------- Equations ---------------------------
% Number of moles of liquid metal (mol)
XM_lSc = (m_lSc/M_Fe) + (m_C/M_C) + (m_Si/M_Si) + (m_Cr/M_Cr) + ...
    (m_Mn/M_Mn) + (m_P / M_P);
% Number of moles in liquid slag zones (mol)
XM_lSl = (m_lSl/M_lSl) + (m_FeO/M_FeO) + (m_SiO2/M_SiO2) + (m_MnO/M_MnO) ...
    + (m_Cr2O3/M_Cr2O3) + (m_P2O5/M_P2O5) + (m_MgO/M_MgO) + (m_CaO/M_CaO);
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

kX_Si = 8.08e-08;
kX_MnO1 = 10;
kX_MnO2 = 10^(2.8*((X_CaO+X_MgO)/X_SiO2)-1.16) * (M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2);
kX_Mn = 189.3;
kX_Cr = 0.3 * ((M_Cr*M_FeO*100) / (M_Cr2O3*M_Fe));
kX_P = 7800;

Xeq_C = 4.9e-4 / X_FeO;
Xeq_Si = kX_Si / (X_FeO^2);
Xeq_MnO1 = X_Mn/kX_MnO1;
% Xeq_MnO2 = sqrt((X_Mn^2*X_SiO2)/(X_Si*kX_MnO2));
Xeq_MnO2 = 0.0001;
Xeq_Mn = X_MnO / (X_FeO * kX_Mn);
Xeq_Cr = X_Cr2O3 / (X_FeO * kX_Cr);
% Xeq_P = sqrt(X_P2O5 / (X_FeO^5 * X_CaO^3 * kX_P));
Xeq_P = 0.0001;

% ------------------- Heat Flows ----------------------
% Energy dissipated from the arcs by conduction (kW)
Q_arc = 0.2 * P_arc;

% Energy flow between the solid and liquid crap zones (kW)
Q_lScsSc = min([m_lSc, m_sSc]) * K_therm1 * K_area1 * (T_lSc - T_sSc);

% Energy flow between solid scrap and solid slag (kW)
Q_sScsSl = min([m_sSc, m_sSl]) * K_therm2 * K_area2 * (T_sSc - T_sSl);

% Energy flow between solid scrap and liquid slag (kW)
Q_sSclSl = min([m_sSc, m_lSl]) * K_therm3 * K_area3 * (T_sSc - T_lSl);

% Energy exchanged between the solid scrap and surrounding gas (kW)
Q_sScgas = (m_sSc/m_EAF) * K_therm4 * (T_sSc-T_gas) * (1-K_sSclSc);

% Portion of solid steel energy lost due to cooling of the furnace walls
Q_sScwater = K_water1 * (T_sSc - T_wall) * (T_sSc/T_melt) * (1 - exp(-(m_sSc/m_EAF)));

% Energy exchange between liquid metal and solid slag (kW)
Q_lScsSl = min([m_lSc, m_sSl]) * K_therm5 * K_area5 * (T_lSc - T_sSl);

% Energy exchange between liquid metal and liquid slag (kW)
Q_lSclSl = min([m_lSc, m_lSl]) * K_therm6 * K_area6 * (T_lSc - T_lSl);

% Energy exchanged between the liquid metal and surrounding gas (kW)
Q_lScgas = (m_lSc/m_EAF) * K_therm7 * (T_lSc - T_gas) * K_sSclSc;

% Energy loss in liquid metal zone due to cooling (kW)
Q_lScwater = K_water2 * (T_lSc - T_wall) * (T_lSc/T_melt) * (1 - exp(-(m_lSc/m_EAF)));

% Energy loss in solid slag due to cooling
Q_sSlwater = K_water3 * (T_sSl - T_wall) * (T_sSl/T_melt) * (1 - exp(-(m_sSl/m_EAF)));

% Net energy in solid slag zone
Q_sSl = Q_sScsSl + Q_lScsSl - Q_sSlwater;

% Heat exchange between liquid slag and gas zone
Q_lSlgas = (m_lSl/m_EAF) * K_therm8 * (T_lSl - T_gas) * K_sSclSc;

% Energy loss in liquid slag due to cooling
Q_lSlwater = K_water4 * (T_lSl - T_wall) * (T_lSl/T_melt) * (1 - exp(-(m_lSl/m_EAF)));

% Net energy in liquid slag zone
Q_lSl = Q_lSclSl + Q_sSclSl - Q_lSlgas - Q_lSlwater;

% Energy received by gas zone from arcs
Q_arcgas = 0.025 * P_arc;

% Energy loss from gas zone to furnace roof and walls
Q_gaswater = K_water5*((T_gas - T_roof)*(A1/(A1+A2)) + (T_gas - T_wall)*(A2/(A1+A2)));

% Radiative heat transfer at arc
Q_arcRAD = 0.75 * P_arc;

% ---------------------- Temperature change -----------------------
% Temperature change of sSl
dT_sSl = (Q_sSl*(1-(T_sSl/T_melt))) / ((m_sSl/M_sSl)*Cp_sSl);
% Temperature change of lSl
dT_lSl = Q_lSl/((m_lSl/M_lSl)*Cp_lSl);

% ------------------------ Mass changes ---------------------------

% Melt rate of solid slag (kg/s)
dm_sSl = -(Q_sSl*(T_sSl/T_melt)) / ((lambda_sSl + Cp_sSl*(T_melt - T_sSl))/M_sSl);
% Liquid slag mass rate (kg/s)
dm_lSl = -dm_sSl;

% Rate of change of injected C
x1_d1 = C_inj;
% x1_d2 by dimensional analysis kg^2/s but in the literature it says kg/s
x1_d2 = -(m_FeO*kd_CL*m_CL)/(m_lSl+m_FeO+m_SiO2+m_MnO+m_Cr2O3+m_P2O5);
% x1_d3 kg but in literature it says kg/s
x1_d3 = -(m_CL*T_lSc*Cp_lSc*(T_air/T_melt)) / (lambda_C + Cp_C*(T_melt-T_air));
dm_CL = x1_d1 + x1_d2 + x1_d3;

% Rate of dissolved C in bath
x2_d1 = -kd_CD*(X_C - Xeq_C); % kg/s
% Should be kg^2/s^2
x2_d2 = -kd_C1*(X_C - Xeq_C) * O2_lance * K_O2CO;
x2_d3 = -x1_d3;
x2_d4 = -kd_Mn1*(M_C/M_MnO)*(X_MnO - Xeq_MnO1); % kg/s
x2_d5 = -kd_C2 * (X_C - Xeq_C) * O2_lance * K_O2CO2;
dm_CD = x2_d1 + x2_d2 + x2_d3 + x2_d4 + x2_d5;

% Rate of change of Si
x3_d1 = -kd_Si1 * (X_Si - Xeq_Si);
x3_d2 = -kd_Si2 * (X_Si - Xeq_Si) * O2_lance * K_O2SiO2;
x3_d3 = -kd_Mn2 * (M_Si/M_MnO) * (X_MnO - Xeq_MnO2);
dm_Si = x3_d1 + x3_d2 + x3_d3;

% Rate of change of Mn
x4_d1 = -kd_Mn * (X_Mn - Xeq_Mn);
x4_d2 = -(M_Mn/M_C) * x2_d4;
x4_d3 = -(M_Mn/M_Si) * x3_d3;
dm_Mn = x4_d1 + x4_d2 + x4_d3;

% Rate of change of Cr
x5_d1 = -2 * kd_Cr1 * (X_Cr - Xeq_Cr);
x5_d2 = -4 * kd_Cr2 * (X_Cr - Xeq_Cr) * O2_lance * K_O2Cr2O3;
dm_Cr = x5_d2 + x5_d1;

% Rate of change of P
dm_P = -2 * kd_P * (X_P - Xeq_P);

% Rate of change of Feo
x7_d1 = (M_FeO/M_C) * x1_d2;
x7_d2 = 2*(M_FeO/M_Si) * x3_d1;
x7_d3 = 2*(M_FeO/M_O2) * O2_lance * K_O2FeO;
x7_d4 = DRI_add * K_FeODRI;
x7_d5 = (M_FeO/M_C) * x2_d1;
x7_d6 = (3/2) * (M_FeO/M_Cr) * x5_d1;
x7_d7 = (M_FeO/M_Mn) * x4_d1;
x7_d8 = (5/2) * (M_FeO/M_P) * dm_P;
dm_FeO = x7_d1 + x7_d2 + x7_d3 + x7_d4 + x7_d5 + x7_d6 + x7_d7 + x7_d8;

% Rate of change of SiO2
dm_SiO2 = -(M_SiO2/M_Si) * dm_Si;

% Rate of change of MnO
dm_MnO = -(M_MnO/M_Mn) * dm_Mn;

% Rate of change of Cr2O3
dm_Cr2O3 = -(M_Cr2O3/(2*M_Cr)) * dm_Cr;

% Rate of change of P2O5
dm_P2O5 = -(M_P2O5/(2*M_P)) * dm_P;

% Rate of change of Fe
x8_d1 = -(M_Fe/M_C) * x1_d2;
x8_d2 = -2 * (M_Fe/M_Si) * x3_d1;
x8_d3 = -(M_Fe/M_FeO) * x7_d3;
x8_d4 = DRI_add * K_FeDRI;
x8_d5 = -(M_Fe/M_C) * x2_d1;
x8_d6 = -(3/2) * (M_Fe/M_Cr) * x5_d1;
x8_d7 = -(M_Fe/M_Mn) * x4_d1;
x8_d8 = -(5/2) * (M_Fe/M_P) * dm_P;
dm_Fe = x8_d1 + x8_d2 + x8_d3 + x8_d4 + x8_d5 + x8_d6 + x8_d7 + x8_d8;

% Rate of change of CO
x9_d1 = -((hd*u1*m_CO)/((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2)));
x9_d2 = -(M_CO/M_C) * (x1_d2 + x2_d1 + x2_d2 + x2_d4);
x9_d3 = 2 * M_CO * k_air1 * k_PR * rp;
x9_d4 = -((k_PR*rp*m_CO)/(m_CO + m_CO2 + m_N2 + m_O2));
x9_d5 = -2*(M_CO/M_O2) * O2_post * K_mCO;
if rp > 0
    dm_CO = x9_d1 + x9_d2 + x9_d4 + x9_d5;
else
    dm_CO = x9_d1 + x9_d2 + x9_d3 + x9_d5;
end

% Rate of change of CO2
x10_d1 = -(hd*u1*m_CO2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x10_d2 = 2*M_CO2*k_air1*k_PR*rp;
x10_d3 = 2*(M_CO2/M_O2)*O2_post*K_mCO;
x10_d4 = (M_CO2/M_CH4)*CH4_inj;
x10_d7 = -(k_PR*rp*m_CO2) / (m_CO+m_CO2+m_N2+m_O2);
x10_d8 = -(M_CO2/M_C) * x2_d5;

% Rate of change of N2
x11_d1 = -(hd*u1*m_N2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x11_d3 = -(k_PR*rp*m_N2) / (m_CO+m_CO2+m_N2+m_O2);
dm_N2 = x11_d1 + x11_d3;

% Rate of change of O2
x12_d1 = -(hd*u1*m_O2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
x12_d4 = O2_post * (1-K_mCO);
x12_d8 = -(M_O2/M_FeO) * x7_d3;


% Electrode Oxidation
dm_el = 3*(R_tip * (I_arc^2/3600) + R_side * (A_side/3600));

% Rate of change of the combustible materials
dm_comb = -kd_comb * m_comb * (T_sSc/T_melt);

% ======================================
% ==== Energy of chemical reactions ====
% ======================================
% All in kW

% a) Fe + 1/2O2 -> FeO
dH_Ta = (x8_d3/M_Fe) * (dH_FeO - dH_FeS + (Cp_FeO - Cp_Fe - 0.5*Cp_O2) * (T_lSc - 298));

% b) FeO + C -> Fe + CO
dH_Tb = (x1_d2+x2_d2)/M_C * (dH_CO - dH_CS - dH_FeO + (Cp_Fe + Cp_CO - Cp_C - Cp_FeO) * (T_lSc - 298));

% c) FeO + Mn -> Fe + MnO
dH_Tc = (x4_d1/M_Mn) * (dH_MnO - dH_FeO - dH_MnS + (Cp_Fe + Cp_MnO - Cp_FeO - Cp_Mn) * (T_lSc - 298));

% d) 2FeO + Si -> 2Fe + SiO2
dH_Td = (x3_d1/M_Si) * ((dH_SiO2+dH_SiO2S-dH_FeO-dH_SiS) + ...
    (2*Cp_Fe + Cp_SiO2 - 2*Cp_FeO - 2*Cp_Si)*(T_lSc-298));

% e) 3FeO + 2Cr -> 3Fe + Cr2O3
dH_Te = (x5_d1/M_Cr) * (dH_Cr2O3 - 3*dH_FeO - dH_CrS + ...
    (3*Cp_Fe + Cp_Cr2O3 - 3*Cp_FeO - 2*Cp_Cr) * (T_lSc-298));

% f) 5FeO + 2P -> 5Fe + P2O5
dH_Tf = (dm_P/M_P) * (dH_P2O5 - 5*dH_FeO - 2*dH_PS + ...
    (5*Cp_Fe + Cp_P2O5 - 5*Cp_FeO - 2*Cp_P) * (T_lSc-298));

% g) C + 1/2O2 -> CO
dH_Tg = (x2_d2/M_C) * ((dH_CO-dH_CS) + (Cp_CO - Cp_C - 0.5*Cp_O2)*(T_lSc-298));

% h) CO + 1/2O2 -> CO2
dH_Th = (x9_d4/M_CO) * ((dH_CO2-dH_CO) + (Cp_CO2 - Cp_CO - 0.5*Cp_O2)*(T_gas-298));

% i) C + O2 -> CO2
dH_Ti = (x2_d5/M_C) * ((dH_CO2-dH_CS) + (Cp_CO2 - Cp_C - Cp_O2)*(T_lSc-298));

% j) MnO + C -> Mn + CO
dH_Tj = (x4_d2/M_Mn) * ((dH_CO+dH_MnS-dH_MnO-dH_CS) + ...
    (Cp_Mn + Cp_CO - Cp_MnO - Cp_C) * (T_lSc - 298));

% k) 2MnO + Si -> 2Mn + SiO2
dH_Tk = (x3_d3/M_Si) * ((dH_SiO2 + dH_SiO2S - dH_MnO - dH_SiS) + ...
    (2*Cp_Mn + Cp_SiO2 - Cp_Si - 2*Cp_MnO) * (T_lSc - 298));

% l) Si + O2 -> SiO2
dH_Tl = (x3_d2/M_Si) * ((dH_SiO2 + dH_SiO2S - dH_SiS) + ...
    (Cp_SiO2 - Cp_Si - 2*Cp_O2) * (T_lSc - 298));

% m) 2Cr + 3/2O2 -> Cr2O3
dH_Tm = (x5_d2/M_Cr) * ((dH_Cr2O3 - 2*dH_CrS) + ...
    (Cp_Cr2O3 - 2*Cp_Cr - 0.5*Cp_O2) * (T_lSc - 298));

% n) CH4 + 2O2 -> CO2 + 2H2O
dH_Tn = -(CH4_inj/M_CH4) * ((dH_CO2 + 2*dH_H2O - dH_CH4) + ...
    (Cp_CO2 + 2*Cp_H2O - Cp_CH4 - 2*Cp_O2) * (T_gas - 298));

% o) Graphite to CO2
dH_To = -(dm_el/M_C) * (dH_CO2 + (Cp_CO2 - Cp_C - Cp_O2) * (T_gas - 298));

% p) C9H20 + 54O2 -> 9CO2 + 90H2O
dH_Tp = -(dm_comb/M_C9H20) * ((dH_CO2 + dH_H2O - dH_C9H20) + ...
    (Cp_CO2 + Cp_H2O - Cp_C9H20 - Cp_O2) * (T_gas - 298));

% Total energy of the chemical reactions for liquid metal
Q_lScchem = dH_Ta + dH_Tb + dH_Tc + dH_Td + dH_Te + dH_Tf + dH_Tg ...
    + dH_Ti + dH_Tj + dH_Tk + dH_Tl + dH_Tm + dH_Tn + dH_To + dH_Tp;

% Energy added to solid steel zone from the oxygen burners
Q_CH4 = dH_Tn * K_burn * (0.35 + 0.65*tanh((1573/T_sSc)-1));

% Energy received from oxygen burners
Q_CH4gas = dH_Tn * (1 - K_burn * (0.35 + 0.65*tanh((1573/T_sSc)-1)));

% ==================== Dependent Equations =====================

x10_d5 = (M_CO2/M_C) * dm_el;
x10_d6 = -9 * (M_CO2/M_C9H20)*dm_comb;
if rp > 0
    dm_CO2 = x10_d8 + x10_d7 + x10_d6 + x10_d5 + x10_d4 + x10_d3 + x10_d1;
else
    dm_CO2 = x10_d8 + x10_d2 + x10_d6 + x10_d5 + x10_d4 + x10_d3 + x10_d1;
end

% Rate of change of O2
x12_d2 = -M_O2 * k_air1 * k_PR * rp;
x12_d3 = -(M_O2/(2*M_CO2)) * x10_d2;
x12_d5 = -(M_O2/M_C) * dm_el;
x12_d6 = -14 * (M_O2/M_C9H20) * dm_comb;
x12_d7 = -(k_PR*rp*m_O2) / (m_CO+m_CO2+m_N2+m_O2);
if rp > 0
    dm_O2 = x12_d8 + x12_d7 + x12_d6 + x12_d5 + x12_d4 + x12_d3 + x12_d1;
else
    dm_O2 = x12_d8 + x12_d2 + x12_d6 + x12_d5 + x12_d4 + x12_d3 + x12_d1;
end

% Gas zone energy balance
Q_gas = Q_arcgas + (1-K_post)*dH_Th + Q_CH4gas + Q_sScgas + Q_lScgas ...
    + Q_lSlgas - Q_gaswater;

% Temperature change of gas
dT_gas = Q_gas/((m_gas/M_gas)*Cp_gas);

% Relative Pressure
dp = (R*T_gas/V_gas) * (dm_CO/M_CO + dm_CO2/M_CO2 + dm_N2/M_N2 + dm_O2/M_O2) * tc ...
    + (R*dT_gas*tc/V_gas) * (m_CO/M_CO + m_CO2/M_CO2 + m_N2/M_N2 + m_O2/M_O2);

% Energy balance for the solid steel zone (sSc)
Q_sSc = (Q_arc + Q_CH4 + K_post*dH_Th)*(1-K_sSclSc) + Q_lScsSc ...
    - Q_sScsSl - Q_sSclSl - Q_sScgas - Q_sScwater - Q_sScRAD;

% Energy received by liquid metal zone (lSc)
Q_lSc = (Q_arc + Q_CH4 + K_post*dH_Th)*K_sSclSc + Q_lScchem ...
    - Q_lScsSc - Q_lScsSl - Q_lSclSl - Q_lScgas - Q_lScwater - Q_lScRAD;

% Melt rate of solid scrap
dm_sSc = -((Q_sSc*(T_sSc/T_melt)) / (lambda_sSc + Cp_sSc*(T_melt - T_sSc))) * M_Fe;

% Liquid metal mass rate
dm_lSc = -dm_sSc;

% Radiosity of roof
J_roof = (ep1*sig*T_roof^4/1000 - (1-ep1)*(VF12*J_wall + VF13*J_sSc ...
    + VF14*J_lSc + VF15*Q_arcRAD));

% Radiosity of wall
J_wall = (ep2*sig*T_wall^4/1000 - (1-ep2)*(VF21*J_wall + VF23*J_sSc ...
    + VF24*J_lSc + VF25*Q_arcRAD));

% Radiosity of sSc
J_sSc = (ep3*sig*T_sSc^4/1000 - (1-ep3)*(VF31*J_roof + VF32*J_wall ...
    + VF35*Q_arcRAD));

% Radiosity of lSc
J_lSc = (ep4*sig*T_lSc^4/1000 - (1-ep4)*(VF41*J_roof + VF42*J_wall ...
    + VF45*Q_arcRAD));

% Radiative heat flow in roof
Q_roofRAD = A1 * (VF12*(J_roof-J_wall) + VF13*(J_roof-J_sSc) ...
    + VF14*(J_roof-J_lSc)) - VF51 * Q_arcRAD;

% Radiative heat flow in wall
Q_wallRAD = A2 * (VF21*(J_wall-J_roof) + VF23*(J_wall-J_sSc) ...
    + VF24*(J_wall-J_lSc)) - VF52 * Q_arcRAD;

% Radiative heat flow in sSc
Q_sScRAD = A3 * (VF31*(J_sSc-J_roof) + VF32*(J_sSc-J_wall)) ...
    - VF53 * Q_arcRAD;

% Radiative heat flow in lSc
Q_lScRAD = A4 * (VF41*(J_lSc-J_roof) + VF42*(J_lSc-J_wall)) ...
    - VF54 * Q_arcRAD;

% Temperature change of sSc
dT_sSc = (Q_sSc*(1-(T_sSc/T_melt))) / ((m_sSc/M_Fe)*Cp_sSc);

% Temperature change of lSc
dT_lSc = Q_lSc/((m_lSc/M_Fe)*Cp_lSc);

% Temperature change of roof
dT_roof = (-Q_roofRAD + (A1/(A1+A2))*Q_gaswater - phi1*Cp_H2O*(T_roof - T_water)) ...
    / (A1*d1*rho*Cp_roof);

% Temperature change of wall
dT_wall = (-Q_wallRAD + (A2/(A1+A2))*Q_gaswater - phi2*Cp_H2O*(T_wall - T_water)) ...
    / (A2*d2*rho*Cp_wall);

% ======================= Update value ============================
m_CL = m_CL + dm_CL * tc;
m_C = m_C + (dm_CL + dm_CD) * tc;
m_CO = m_CO + dm_CO * tc;
m_CO2 = m_CO2 + dm_CO2 * tc;
m_comb = m_comb + dm_comb * tc;
m_Cr = m_Cr + dm_Cr * tc;
m_Cr2O3 = m_Cr2O3 + dm_Cr2O3;
m_Fe = m_Fe + dm_Fe * tc;
m_FeO = m_FeO + dm_FeO * tc;
m_lSc = m_lSc + dm_lSc * tc;
m_lSl = m_lSl + dm_lSl * tc;
m_Mn = m_Mn + dm_Mn * tc;
m_MnO = m_MnO + dm_MnO * tc;
m_N2 = m_N2 + dm_N2 * tc;
m_O2 = m_O2 + dm_O2 * tc;
m_P = m_P + dm_P * tc;
m_P2O5 = m_P2O5 + dm_P2O5 * tc;
m_sSc = m_sSc + dm_sSc * tc;
m_sSl = m_sSl + dm_sSl * tc;
dm_gas = dm_CO + dm_CO2 + dm_N2 + dm_O2;
m_gas = m_gas + dm_gas * tc;
m_Si = m_Si + dm_Si * tc;
m_SiO2 = m_SiO2 + dm_SiO2 * tc;

T_gas = T_gas + dT_gas * tc;
T_lSc = T_lSc + dT_lSc * tc;
T_lSl = T_lSl + dT_lSl * tc;
T_roof = T_roof + dT_roof * tc;
T_sSc = T_sSc + dT_sSc * tc;
T_sSl = T_sSl + dT_sSl * tc;
T_wall = T_wall + dT_wall * tc;

% rp = rp + dp * tc;

if mod(step, 1/tc) == 0
    gas_temp(step*tc) = T_gas;
    sSc_temp(step*tc) = T_sSc;
    sSl_temp(step*tc) = T_sSl;
    lSc_temp(step*tc) = T_lSc;
    lSl_temp(step*tc) = T_lSl;
    steel_Fe(step*tc) = X_lSc;
end

if T_gas > 2000
    CH4_inj = 0;
else
    CH4_inj = 0.1;
end

end

% Graph generation
% time = linspace(1, secs, secs);
% 
% figure
% plot(time, gas_temp)
% hold on
% plot(time, sSc_temp)
% plot(time, sSl_temp)
% plot(time, lSc_temp)
% plot(time, lSl_temp)
% 
% legend('Gas', 'Solid Metal', 'Solid Slag', 'Liquid Metal', 'Liquid Slag')
% hold off
% 
% figure
% plot(time, steel_Fe)