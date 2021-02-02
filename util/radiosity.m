% Initial radiosity assumption
J_roof = -20000;
J_wall = -7700;
J_sSc = 110000;
J_lSc = 240000;

% Radiosity of roof
J_roof = ep1*sig*T_roof^4 - (1-ep1)*(VF12*J_wall + VF13*J_sSc ...
    + VF14*J_lSc + VF15*Q_arcRAD);

% Radiosity of wall
J_wall = ep2*sig*T_wall^4 - (1-ep2)*(VF21*J_wall + VF23*J_sSc ...
    + VF24*J_lSc + VF25*Q_arcRAD);

% Radiosity of sSc
J_sSc = ep3*sig*T_sSc^4 - (1-ep3)*(VF31*J_roof + VF32*J_wall ...
    + VF35*Q_arcRAD);

% Radiosity of lSc
J_lSc = ep4*sig*T_lSc^4 - (1-ep4)*(VF41*J_roof + VF42*J_wall ...
    + VF45*Q_arcRAD);

J_roof
J_wall
J_sSc
J_lSc