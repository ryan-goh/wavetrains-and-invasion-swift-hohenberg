function [dydt] = tw_vec(y,c,m)
    u = y(1);
    v = y(2);
    dudt = v
    dvdt = -(m*u - u.^2 +c*v);

    dydt = [dudt;dvdt];
end