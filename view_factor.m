% -------------------------- Geometry ------------------------------
% 1: Roof = Disk
% 2: Wall = Cylinder
% 3: steel scrap = Disk
% 4: Liquid Metal = Disk
% 5: Arc = Cylinder

clc
clear

r_eafout = 6.6/2; % m
r_hole = 3.4/2; % m
r_eafin = 4.9/2;
h_eafup = 2.9;
h_eaflow = 1;
h_wall = 0.8247;
h_arc = 0.6; % Arbitrary

h_bath = h_eaflow/2;
h_scrap = h_eafup + h_eaflow - h_wall;
h_cone = h_arc;

K_sSclSc = 0.5 * tanh(5*(h_bath-h_scrap+h_cone)) + 0.5;

A1 = (r_eafout^2 - r_hole^2) * pi;
A2 = pi*r_eafout*2*h_wall;
A3 = (r_eafout^2 - r_hole^2) * pi;
A4 = pi*r_eafin^2;
A5 = pi*r_hole*2*h_arc;

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
VF43 = (A3/A4) * VF34;
VF45 = (A5/A4) * VF54;
