
% FeO content in DRI
DRI_FeO = 0.07; 
% ------------------------ Thermal Constants -------------------------
% Heat transfer coeff. (can be changed depending on the observed melting
% time)
k_ther1 = 0.4; % kW.K-1.m-2
k_ther5 = 0.2;

% -------------------------- Air Constants ---------------------------
k_air1 = 7.3; % mol.kg-1
k_air2 = 27.4; % mol.kg-1

% ------------------------- Area Constants ---------------------------
% Specific area (0.005m^2.kg-1 corresponds to sheets of thickness 0.05m)
k_area1 = 0.005; % m^2.kg-1
% Correspeonds to cubic blocks of density 2500kg.m-3 with 0.02m sides
k_area2 = 0.12;

% ------------------------- Rate Constants ---------------------------
k_dC = 12; % kg.s-1
k_dSi = 144; % kg.s-1

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
% Molecular weight of Slag


% --------------------- Dimensionless Constants -----------------------
k_gr = 15.0;
k_u = 8.44;
k_xc = 491e-06;
k_xsi = 8.08e-08;

% --------------------- Thermodynamic Properties -----------------------
% Latent heat of fusion of Fe
lambda_Fe = 13794; % kJ.mol-1
% Specific heat capacity of Fe (Currently at room temp.)
% TODO: Change the heat capacity to operating temperature
Cp_Fes = 24.795; % kJ.mol-1.K-1

% ------------------------------ ODEs ----------------------------------
syms x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t)

% Rate of Change of Solid Scrap (kg.s-1)
ode1 = diff(x1) == (-M_Fe*k_ther1*k_area1*x1*(x12-x13)*(x13/x12))/((lambda_Fe+Cp_Fes)*(x12-x13));
% Rate of Change of Liquid Metal (kg.s-1)
ode2 = diff(x2) == (M_Fe*k_ther1*k_area1*x1*(x12-x13)*(x13/x12))/((lambda_Fe+Cp_Fes)*(x12-x13)) ...
    + (x7*k_gr*M_Fe*d5)/((x6+x7+x8)*M_C) + (M_Fe/M_C)*k_dC*(X_C-Xeq_C) ...
    + (2*M_Fe/M_Si)*k_dSi*(X_Si - Xeq_Si) - (2*M_Fe*d1/M_O2) + 0.93*d2;
% Rate of change of dissolved carbon in steel melt
ode3 = diff(x3) == -k_dC*(X_C - Xeq_C);
% Rate of decrease of silicon in the steel melt
ode4 = diff(x4) == -k_dSi*(X_Si - Xeq_Si);
% Rate of change of solid slag
ode5 = diff(x5) == (-M_slag*k_ther5*k_area5*x5*(x12-x13)*(x13/x12))/((lambda_slag+Cp_slag)*(x12-x13)) + d3;
% Rate of change of liquid slag
ode6 = diff(x6) == (M_slag*k_ther5*k_area5*x5*(x12-x13)*(x13/x12))/((lambda_slag+Cp_slag)*(x12-x13));
% Rate of change of FeO in slag
ode7 = diff(x7) == (2*M_FeO*d1/M_O2) - (x7*k_gr*M_FeO*d5)/((x6+x7+x8)*M_C) ...
    + DRI_FeO*d2 - (M_FeO/M_C)*k_dC*(X_C - Xeq_C) - (2*M_FeO/M_Si)*k_dSi*(X_Si - Xeq_Si);
% Rate of change of SiO2 in slag
ode8 = diff(x8) == (M_SiO2/M_Si)*k_dSi*(X_Si - Xeq_Si);
% Rate of change of CO in gas phase
% Assumed positive relative pressure
ode9 = diff(x9) == (M_CO*d5/M_C) - (h_d*u1*x9)/((k_u*u2+h_d)*(x9+x10+x11)) ...
    + (M_CO/M_C)*k_dC*(X_C - Xeq_C) + (k_PR*x14*x9)/(x9+x10+x11);
% Rate of change of CO2 in gas phase
ode10 = diff(x10) == -(h_d*u1*x10)/((k_u*u2+h_d)*(x9+x10+x11)) ...
    - (k_PR*x14*x10)/(x9+x10+x11);
% Rate of change of N2 in gas phase
ode11 = diff(x11) == -(h_d*u1*x11)/((k_u*u2+h_d)*(x9+x10+x11)) ...
    - (k_PR*x14*x11)/(x9+x10+x11) + d5/150;
% Temperature of liquid metal
ode12 = diff(x12) == (pt + d4 - k_VT*

% ---------------------------- Equations -------------------------------
X_C = (x3/M_C) / ((x2/M_Fe) + (x3/M_C) + (x4/M_Si));
Xeq_C = k_xc*((x6*M_FeO)/(x7*M_slag) + (x8*M_FeO)/(x7*M_SiO2) + 1);
X_Si = (x4/M_Si) / ((x2/M_Fe) + (x3/M_C) + (x4/M_Si));
Xeq_Si = k_xsi*((x6*M_FeO)/(x7*M_slag) + (x8*M_FeO)/(x7*M_SiO2) + 1)^2;

