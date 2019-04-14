function res = Sigma( x, u, s )
    sigma0 = 1.63e-3;
    t_crit = 1414;
    teta = u + t_crit;
    if (x < s) %hot [liquid]
        res = sigma0 * exp(u/teta);
    else %cold [solid]
        res = teta^2 / (teta^3 + 1);
    end;
end

