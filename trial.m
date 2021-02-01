fun = @eaf_improved;
x0 = linspace(1, 28, 28);
options = optimset('disp', 'iter', 'LargeScale', 'off', 'TolFun', .001 ...
    , 'MaxIter', 10000000, 'MaxFunEvals', 100000);
dx = fsolve(fun, x0, options);