fun = @eaf_improved;
% Initial value
x0 = [
    % Q_sSc
    50000,
    % Q_lSc
    0,
    % J_roof
    0,
    % J_wall
    0,
    % J_sSc
    0,
    % J_lSc
    0,
    % Q_roofRAD
    -5000,
    % Q_wallRAD
    -6000,
    % Q_sScRAD
    -11000,
    % Q_lScRAD
    0,
    % dT_sSc
    4,
    % dT_lSc
    3,
    % dT_roof
    0.5,
    % dT_wall
    0.5,
    % dm_sSc
    -500,
    % dm_lSc
    500
    
    ];
options = optimset('Algorithm', 'levenberg-marquardt', 'disp', 'iter', 'LargeScale', 'off', 'TolFun', .001 ...
    , 'MaxIter', 10000000, 'MaxFunEvals', 100000);
dx = fsolve(fun, x0, options);