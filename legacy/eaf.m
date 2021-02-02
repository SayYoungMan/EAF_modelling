function F = eaf(x)
%% Pre-calculation Equations and variables
% ---------------------- Initial Conditions -----------------------
% Arc Current
I_arc = 44; % kA

% Initial Mass (kg)
m_lSc = 106000;
m_lSl = 1000;
m_sSc = 40000;
m_sSl = 8000;
m_gas = 100;

m_C = 160;
m_CL = 0;
m_Si = 240;
m_Cr = 80;
m_Mn = 240;
m_P = 20;
m_FeO = 20;
m_SiO2 = 63;
m_MnO = 10;
m_Cr2O3 = 10;
m_P2O5 = 10;
m_CaO = 5103;
m_MgO = 3708;
m_CO = 0;
m_CO2 = 1;
m_N2 = 10;
m_O2 = 5;

m_comb = 440;

% C Injection rate kg/s
C_inj = 10;

% Oxygen burner rate kg/s
CH4_inj = 0.5;

% O2 lance rate kg/s
O2_lance = 1;

% O2 rate for CO post-combustion kg/s
O2_post = 0;
K_mCO = 0;

% DRI Addition rate kg/s
DRI_add = 1000;

% Initial Temperature
T_sSc = 1245;
T_lSc = 1800;
T_sSl = 1245;
T_lSl = 1800;
T_gas = 800;
T_wall = 300;
T_roof = 300;

% ------------------------- Variables -----------------------------
% EAF mass capacity
m_EAF = 105; % kg
K_post = 100;

% Arc Power kW
P_arc = 3000;

% ------------------------- Constants -----------------------------
% Dimensionless constant for approximation
k_U = 6.44;

% Ratio between reaction rate and relative pressure
k_PR = 0.6;

% Universal Gas Constant
R = 8.314; % J/mol K

% Gas zone volume
V_gas = 45; % m^3

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

M_lSl = 0.05;

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
hd = 0.65;
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

VF51 = (A1/A5) * VF15;
VF52 = (A2/A5) * VF25;
VF53 = (1-VF51-VF52) * K_sSclSc;
VF54 = (1-VF51-VF52) * (1-K_sSclSc);

VF31 = VF13;
VF32 = VF12;
VF34 = (1/2)*(S3 - sqrt(S3^2 - 4*(r_eafout/ratio/r_eafin)^2))/ratio;
VF35 = 1 - VF31 - VF32 - VF34;

VF41 = (A1/A4) * VF14;
% VF43 = (A3/A4) * VF34;
VF45 = (A5/A4) * VF54;

% ----------------------------- Equations ---------------------------
% Number of moles of liquid metal
XM_lSc = (m_lSc/M_Fe) + (m_C/M_C) + (m_Si/M_Si) + (m_Cr/M_Cr) + ...
    (m_Mn/M_Mn) + (m_P / M_P);
% Number of moles in liquid slag zones
XM_lSl = (m_lSl/M_lSl) + (m_FeO/M_FeO) + (m_SiO2/M_SiO2) + (m_MnO/M_MnO) ...
    + (m_Cr2O3/M_Cr2O3) + (m_P2O5/M_P2O5);
% Mole fraction of FeO
X_C = (m_C/M_C) / XM_lSc;
X_Si = (m_Si/M_Si) / XM_lSc;
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

kX_Si = 5.7e-4 * ((M_lSl^2*M_Fe) / (M_FeO^2 * M_Si * 100^3));
kX_MnO1 = 0.064 * X_MnO;
kX_MnO2 = 10^(2.8*((X_CaO+X_MgO)/X_SiO2)-1.16) * (M_MnO^2 * M_Si * M_Fe) / (M_Mn^2 * M_lSl * M_SiO2);
kX_Mn = 1.9 * ((M_FeO * M_Mn * 100) / (M_MnO*M_Fe));
kX_Cr = 0.3 * ((M_Cr*M_FeO*100) / (M_Cr2O3*M_Fe));
kX_P = 9.8e-4 * ((M_Cr * M_FeO * 100) / (M_Cr2O3 * M_Fe));

Xeq_C = 4.9e-4 / X_FeO;
Xeq_Si = kX_Si / (X_FeO^2);
Xeq_MnO1 = X_Mn/kX_MnO1;
Xeq_MnO2 = sqrt((X_Mn^2*X_SiO2)/(X_Si*kX_MnO2));
Xeq_Mn = X_MnO / (X_FeO * kX_Mn);
Xeq_Cr = X_Cr2O3 / (X_FeO * kX_Cr);
Xeq_P = sqrt(X_P2O5 / (X_FeO^5 * X_CaO^3 * kX_P));
%% System of equations
% ------------------------- System of equations ---------------------

% List of variables
% x(1) = Q_sSc
% x(2) = Q_arc
% x(3) = Q_CH4
% x(4) = Q_lScsSc
% x(5) = Q_sScsSl
% x(6) = Q_sSclSl
% x(7) = Q_sScgas
% x(8) = Q_sScwater
% x(9) = Q_lSc
% x(10) = Q_lScsSl
% x(11) = Q_lSclSl
% x(12) = Q_lScgas
% x(13) = Q_lScwater
% x(14) = Q_sSl
% x(15) = Q_sSlwater
% x(16) = Q_lSl
% x(17) = Q_lSlgas
% x(18) = Q_lSlwater
% x(19) = Q_gas
% x(20) = Q_arcgas
% x(21) = Q_CH4gas
% x(22) = Q_gaswater
% x(23) = J_roof = J1
% x(24) = J_wall = J2
% x(25) = J_sSc = J3
% x(26) = J_lSc = J4
% x(27) = J_arc = J5 = Q_arcRAD
% x(28) = Q_roofRAD
% x(29) = Q_wallRAD
% x(30) = Q_sScRAD
% x(31) = Q_lScRAD
% x(32) = dT_sSc
% x(33) = dT_lSc
% x(34) = dT_sSl
% x(35) = dT_lSl
% x(36) = dT_gas
% x(37) = dT_roof
% x(38) = dT_wall
% x(39) = dm_sSc
% x(40) = dm_lSc
% x(41) = dm_sSl
% x(42) = dm_lSl

% Energy balance for the solid steel zone (sSc)
F(1) = -x(1) + (x(2) + x(3) + K_post*x(121))*(1-K_sSclSc) + x(4) ...
    - x(5) - x(6) - x(7) - x(8) - x(30);

% Energy dissipated from the arcs by conduction
F(2) = x(2) - 0.2*P_arc;

% Energy added to solid steel zone from the oxygen burners
F(3) = x(3) - x(127) * K_burn * (0.35 + 0.65*tanh((1573/T_sSc)-1));

% Energy flow between the solid and liquid crap zones
F(4) = x(4) - min([m_lSc, m_sSc]) * K_therm1 * K_area1 * (T_lSc - T_sSc);

% Energy flow between solid scrap and solid slag
F(5) = x(5) - min([m_sSc, m_sSl]) * K_therm2 * K_area2 * (T_sSc - T_sSl);

% Energy flow between solid scrap and liquid slag
F(6) = x(6) - min([m_sSc, m_lSl]) * K_therm3 * K_area3 * (T_sSc - T_lSl);

% Energy exchanged between the solid scrap and surrounding gas
F(7) = x(7) - (m_sSc/m_EAF) * K_therm4 * (T_sSc-T_gas) * (1-K_sSclSc);

% Portion of solid steel energy lost due to cooling of the furnace walls
F(8) = x(8) - K_water1 * (T_sSc - T_wall) * (T_sSc/T_melt) ...
    * (1 - exp(-(m_sSc/m_EAF)));

% Energy received by liquid metal zone (lSc)
F(9) = -x(9) + (x(2) + x(3) + K_post*x(121))*K_sSclSc + x(130) ...
    - x(4) - x(10) - x(11) - x(12) - x(13) - x(31);

% Energy exchange between liquid metal and solid slag
F(10) = x(10) - min([m_lSc, m_sSl]) * K_therm5 * K_area5 * (T_lSc - T_sSl);

% Energy exchange between liquid metal and liquid slag
F(11) = x(11) - min([m_lSc, m_lSl]) * K_therm6 * K_area6 * (T_lSc - T_lSl);

% Energy exchanged between the liquid metal and surrounding gas
F(12) = x(12) - (m_lSc/m_EAF) * K_therm7 * (T_lSc - T_gas) * K_sSclSc;

% Energy loss in liquid metal zone due to cooling
F(13) = x(13) - K_water2 * (T_lSc - T_wall) * (T_lSc/T_melt) ...
    * (1 - exp(-(m_lSc/m_EAF)));

% Net energy in solid slag zone
F(14) = -x(14) + x(5) + x(10) - x(15);

% Energy loss in solid slag due to cooling
F(15) = x(15) - K_water3 * (T_sSl - T_wall) * (T_sSl/T_melt) ...
    * (1 - exp(-(m_sSl/m_EAF)));

% Net energy in liquid slag zone
F(16) = -x(16) + x(11) + x(6) - x(17) - x(18);

% Heat exchange between liquid slag and gas zone
F(17) = (m_lSl/m_EAF) * K_therm8 * (T_lSl - T_gas) * K_sSclSc;

% Energy loss in liquid slag due to cooling
F(18) = x(18) - K_water4 * (T_lSl - T_wall) * (T_lSl/T_melt) ...
    * (1 - exp(-(m_lSl/m_EAF)));

% Gas zone energy balance
F(19) = -x(19) + x(20) + (1-K_post)*x(121) + x(21) + x(7) + x(12) ...
    + x(17) - x(22);

% Energy received by gas zone from arcs
F(20) = x(20) - 0.025 * P_arc;

% Energy received from oxygen burners
F(21) = x(21) - x(127) * (1 - K_burn * (0.35 + 0.65*tanh((1573/T_sSc)-1)));

% Energy loss from gas zone to furnace roof and walls
F(22) = x(22) -K_water5*((T_gas - T_roof)*(A1/(A1+A2)) + (T_gas - T_wall)*(A2/(A1+A2)));

% Radiative heat transfer from roof
F(23) = x(23) - ep1*sig*T_roof^4 - (1-ep1)*(VF12*x(24) + VF13*x(25) ...
    + VF14*x(26) + VF15*x(27));

% Radiative heat transfer from wall
F(24) = x(24) - ep2*sig*T_wall^4 - (1-ep2)*(VF21*x(23) + VF23*x(25) ...
    + VF24*x(26) + VF25*x(27));

% Radiative heat transfer from Solid Metal
F(25) = x(25) - ep3*sig*T_sSc^4 - (1-ep3)*(VF31*x(23) + VF32*x(24) ...
    + VF35*(27));

% Radiative heat transfer from liquid metal
F(26) = x(26) - ep4*sig*T_lSc^4 - (1-ep4)*(VF41*x(23) + VF42*x(24) ...
    + VF45*(27));

% Radiative heat transfer at arc
F(27) = x(27) - 0.75 * P_arc;

% Radiative heat flow in roof
F(28) = -x(28) + A1 * (VF12*(x(23)-x(24)) + VF13*(x(23)-x(25)) ...
    + VF14*(x(23)-x(26))) - VF51 * x(27);

% Radiative heat flow in wall
F(29) = -x(29) + A2 * (VF21*(x(24)-x(23)) + VF23*(x(24)-x(25)) ...
    + VF24*(x(24)-x(26))) - VF52 * x(27);

% Radiative heat flow in sSc
F(30) = -x(30) + A3 * (VF31*(x(25)-x(23)) + VF32*(x(25)-x(24))) ...
    - VF53 * x(27);

% Radiative heat flow in sSc
F(31) = -x(31) + A4 * (VF41*(x(26)-x(23)) + VF42*(x(26)-x(24))) ...
    - VF54 * x(27);

% Temperature change of sSc
F(32) = x(32) - (x(1)*(1-(T_sSc/T_melt))) / (m_sSc*Cp_sSc);

% Temperature change of lSc
F(33) = x(33) - x(9)/(m_lSc*Cp_lSc);

% Temperature change of sSl
F(34) = x(34) - (x(14)*(1-(T_sSl/T_melt))) / (m_sSl*Cp_sSl);

% Temperature change of lSl
F(35) = x(35) - x(16)/(m_lSl*Cp_lSl);

% Temperature change of gas
F(36) = x(36) - x(19)/(m_gas*Cp_gas);

% Temperature change of roof
F(37) = -x(37) + (-x(28) + (A1/(A1+A2))*x(22) - phi1*Cp_H2O*(T_roof - T_water)) ...
    / (A1*d1*rho*Cp_roof);

% Temperature change of wall
F(38) = -x(38) + (-x(29) + (A2/(A1+A2))*x(22) - phi2*Cp_H2O*(T_wall - T_water)) ...
    / (A2*d2*rho*Cp_wall);

% Melt rate of solid scrap
F(39) = x(39) + (x(1)*(T_sSc/T_melt)) / (lambda_sSc + Cp_sSc*(T_melt - T_sSc));

% Liquid metal mass rate
F(40) = x(40) + x(39);

% Melt rate of solid slag
F(41) = x(41) + (x(14)*(T_sSl/T_melt)) / (lambda_sSl + Cp_sSl*(T_melt - T_sSl));

% Liquid slag mass rate
F(42) = x(42) + x(41);

% Rate of change of injected C
F(43) = x(43) - C_inj;
F(44) = x(44) + (m_FeO*kd_CL*m_CL)/(m_lSl+m_FeO+m_SiO2+m_MnO+m_Cr2O3+m_P2O5);
F(45) = x(45) + (m_CL*T_lSc*Cp_lSc*(T_air/T_melt)) / (lambda_C + Cp_C*(T_melt-T_air));
F(46) = x(46) - x(43) - x(44) - x(45);

% Rate of dissolved C in bath
F(47) = x(47) + kd_CD*(X_C - Xeq_C);
F(48) = x(48) + kd_C1*(X_C - Xeq_C) * O2_lance * K_O2CO;
F(49) = x(49) + x(45);
F(50) = x(50) + kd_Mn1*(M_C/M_MnO)*(X_MnO - Xeq_MnO1);
F(51) = x(51) + kd_C2 * (X_C - Xeq_C) * O2_lance * K_O2CO2;
F(52) = x(52) - x(47) - x(48) - x(49) - x(50) - x(51) - x(52);

% Rate of change of Si
F(53) = x(53) + kd_Si1 * (X_Si - Xeq_Si);
F(54) = x(54) + kd_Si2 * (X_Si - Xeq_Si) * O2_lance * K_O2SiO2;
F(55) = x(55) + kd_Mn2 * (M_Si/M_MnO) * (X_MnO - Xeq_MnO2);
F(56) = x(56) - x(53) - x(54) - x(55);

% Rate of change of Mn
F(57) = x(57) + kd_Mn * (X_Mn - Xeq_Mn);
F(58) = x(58) + (M_Mn/M_C) * x(50);
F(59) = x(59) + (M_Mn/M_Si) * x(55);
F(60) = x(60) - x(57) - x(58) - x(59);

% Rate of change of Cr
F(61) = x(61) + 2*kd_Cr1 * (X_Cr - Xeq_Cr);
F(62) = x(62) + 4*kd_Cr2 * (X_Cr - Xeq_Cr) * O2_lance * K_O2Cr2O3;
F(63) = x(63) - x(62) - x(61);

% Rate of change of P
F(64) = x(64) + 2*kd_P * (X_P - Xeq_P);

% Rate of change of Feo
F(65) = x(65) - (M_FeO/M_C) * x(44);
F(66) = x(66) - 2*(M_FeO/M_Si) * x(53);
F(67) = x(67) - 2*(M_FeO/M_O2) * O2_lance * K_O2FeO;
F(68) = x(68) - DRI_add * K_FeODRI;
F(69) = x(69) - (M_FeO/M_C) * x(47);
F(70) = x(70) - (3/2) * (M_FeO/M_Cr) * x(61);
F(71) = x(71) - (M_FeO/M_Mn) * x(57);
F(72) = x(72) - (5/2) * (M_FeO/M_P) * x(64);
F(73) = x(73) - x(72) - x(71) - x(70) - x(69) - x(68) - x(67) - x(66) - x(65);

% Rate of change of SiO2
F(74) = x(74) + (M_SiO2/M_Si) * x(56);

% Rate of change of MnO
F(75) = x(75) + (M_MnO/M_Mn) * x(60);

% Rate of change of Cr2O3
F(76) = x(76) + (M_Cr2O3/(2*M_Cr)) * x(63);

% Rate of change of P2O5
F(77) = x(77) + (M_P2O5/(2*M_P)) * x(64);

% Rate of change of Fe
F(78) = x(78) + (M_Fe/M_C) * x(44);
F(79) = x(79) + 2*(M_Fe/M_Si) * x(53);
F(80) = x(80) + (M_Fe/M_FeO) * x(67);
F(81) = x(81) - DRI_add * K_FeDRI;
F(82) = x(82) + (M_Fe/M_C) * x(47);
F(83) = x(83) + (3/2) * (M_Fe/M_Cr) * x(61);
F(84) = x(84) + (M_Fe/M_Mn) * x(57);
F(85) = x(85) + (5/2) * (M_Fe/M_P) * x(64);
F(86) = x(86) - x(85) - x(84) - x(83) - x(82) - x(81) - x(80) - x(79) - x(78);

% Rate of change of CO
F(87) = x(87) + ((hd*u1*m_CO)/((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2)));
F(88) = x(88) + (M_CO/M_C) * (x(44) + x(47) + x(48) + x(50));
F(89) = x(89) + ((k_PR*x(111)*m_CO)/(m_CO + m_CO2 + m_N2 + m_O2));
F(90) = x(90) + 2*(M_CO/M_O2) * O2_post * K_mCO;
F(91) = x(91) - x(90) - x(89) - x(88) - x(87);

% Rate of change of CO2
F(92) = x(92) + (hd*u1*m_CO2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
F(93) = x(93) - 2*(M_CO2/M_O2)*O2_post*K_mCO;
F(94) = x(94) - (M_CO2/M_CH4)*CH4_inj;
F(95) = x(95) - (M_CO2/M_C) * x(112);
F(96) = x(96) + 9*(M_CO2/M_C9H20)*x(113);
F(97) = x(97) + (k_PR*x(111)*m_CO2) / (m_CO+m_CO2+m_N2+m_O2);
F(98) = x(98) + (M_CO2/M_C) * x(51);
F(99) = x(99) - x(98) - x(97) - x(96) - x(95) - x(94) - x(93) - x(92);

% Rate of change of N2
F(100) = x(100) + (hd*u1*m_CO2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
F(101) = x(101) + (k_PR*x(111)*m_N2) / (m_CO+m_CO2+m_N2+m_O2);
F(102) = x(102) - x(101) - x(100);

% Rate of change of O2
F(103) = x(103) + (hd*u1*m_O2) / ((k_U*u2+hd)*(m_CO + m_CO2 + m_N2 + m_O2));
F(104) = x(104) + M_O2 * k_air1 * k_PR * x(111);
F(105) = x(105) - O2_post * (1-K_mCO);
F(106) = x(106) + (M_O2/M_C) * x(112);
F(107) = x(107) + 14 * (M_O2/M_C9H20) * x(113);
F(108) = x(108) + (k_PR*x(111)*m_O2) / (m_CO+m_CO2+m_N2+m_O2);
F(109) = x(109) + (M_O2/M_FeO) * x(67);
F(110) = x(110) - x(108) - x(107) - x(106) - x(105) - x(104) - x(103) - x(102);

% Relative Pressure
F(111) = -x(111) + (R*T_gas/V_gas) * (x(91)/M_CO + x(99)/M_CO2 + ...
    x(102)/M_N2 + x(110)/M_O2) + (R*x(36)/V_gas) * (m_CO/M_CO + m_CO2/M_CO2 ...
    + m_N2/M_N2 + m_O2/M_O2);

% Electrode Oxidation
F(112) = x(112) - 3*(R_tip * (I_arc^2/3600) + R_side * (A_side/3600));

% Rate of change of the combustible materials
F(113) = x(113) + kd_comb * m_comb * (T_sSc/T_melt);

% ======================================
% ==== Energy of chemical reactions ====
% ======================================

% a) Fe + 1/2O2 -> FeO
F(114) = x(114) - x(80) * (dH_FeO - dH_FeS + (Cp_FeO - Cp_Fe - 0.5*Cp_O2) ...
    * (T_lSc - 298));

% b) FeO + C -> Fe + CO
F(115) = x(115) - (x(44)+x(48))/M_C * (dH_CO - dH_CS - dH_FeO + ...
    (Cp_Fe + Cp_CO - Cp_C - Cp_FeO) * (T_lSc - 298));

% c) FeO + Mn -> Fe + MnO
F(116) = x(116) - (x(57)/M_Mn) * (dH_MnO - dH_FeO - dH_MnS + (Cp_Fe + Cp_MnO ...
    - Cp_FeO - Cp_Mn) * (T_lSc - 298));

% d) 2FeO + Si -> 2Fe + SiO2
F(117) = x(117) - (x(53)/M_Si) * ((dH_SiO2+dH_SiO2S-dH_FeO-dH_SiS) + ...
    (2*Cp_Fe + Cp_SiO2 - 2*Cp_FeO - 2*Cp_Si)*(T_lSc-298));

% e) 3FeO + 2Cr -> 3Fe + Cr2O3
F(118) = x(118) - (x(61)/M_Cr) * (dH_Cr2O3 - 3*dH_FeO - dH_CrS + ...
    (3*Cp_Fe + Cp_Cr2O3 - 3*Cp_FeO - 2*Cp_Cr) * (T_lSc-298));

% f) 5FeO + 2P -> 5Fe + P2O5
F(119) = x(119) - (x(64)/M_P) * (dH_P2O5 - 5*dH_FeO - 2*dH_PS + ...
    (5*Cp_Fe + Cp_P2O5 - 5*Cp_FeO - 2*Cp_P) * (T_lSc-298));

% g) C + 1/2O2 -> CO
F(120) = x(120) - (x(48)/M_C) * ((dH_CO-dH_CS) + (Cp_CO - Cp_C - 0.5*Cp_O2)*(T_lSc-298));

% h) CO + 1/2O2 -> CO2
F(121) = x(121) - (x(89)/M_CO) * ((dH_CO2-dH_CO) + (Cp_CO2 - Cp_CO - 0.5*Cp_O2)*(T_gas-298));
    
% i) C + O2 -> CO2
F(122) = x(122) - (x(51)/M_C) * ((dH_CO2-dH_CS) + (Cp_CO2 - Cp_C - Cp_O2)*(T_lSc-298));

% j) MnO + C -> Mn + CO
F(123) = x(123) - (x(58)/M_Mn) * ((dH_CO+dH_MnS-dH_MnO-dH_CS) + ...
    (Cp_Mn + Cp_CO - Cp_MnO - Cp_C) * (T_lSc - 298));

% k) 2MnO + Si -> 2Mn + SiO2
F(124) = x(124) - (x(55)/M_Si) * ((dH_SiO2 + dH_SiO2S - dH_MnO - dH_SiS) + ...
    (2*Cp_Mn + Cp_SiO2 - Cp_Si - 2*Cp_MnO) * (T_lSc - 298));

% l) Si + O2 -> SiO2
F(125) = x(125) - (x(54)/M_Si) * ((dH_SiO2 + dH_SiO2S - dH_SiS) + ...
    (Cp_SiO2 - Cp_Si - 2*Cp_O2) * (T_lSc - 298));

% m) 2Cr + 3/2O2 -> Cr2O3
F(126) = x(126) - (x(62)/M_Cr) * ((dH_Cr2O3 - 2*dH_CrS) + ...
    (Cp_Cr2O3 - 2*Cp_Cr - 0.5*Cp_O2) * (T_lSc - 298));

% n) CH4 + 2O2 -> CO2 + 2H2O
F(127) = x(127) + (CH4_inj/M_CH4) * ((dH_CO2 + 2*dH_H2O - dH_CH4) + ...
    (Cp_CO2 + 2*Cp_H2O - Cp_CH4 - 2*Cp_O2) * (T_gas - 298));

% o) Graphite to CO2
F(128) = x(128) + (x(112)/M_C) * (dH_CO2 + (Cp_CO2 - Cp_C - Cp_O2) * (T_gas - 298));

% p) C9H20 + 54O2 -> 9CO2 + 90H2O
F(129) = x(129) + (x(113)/M_C9H20) * ((dH_CO2 + dH_H2O - dH_C9H20) + ...
    (Cp_CO2 + Cp_H2O - Cp_C9H20 - Cp_O2) * (T_gas - 298));

% Total energy of the chemical reactions for liquid metal
F(130) = x(130) - x(114) - x(115) - x(116) - x(117) - x(118) - x(119) - x(120) ...
    - x(122) - x(123) - x(124) - x(125) - x(126);

end

