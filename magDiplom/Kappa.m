function res = Kappa( x, s )
    k1 = 0.56; %28.21
    k2 = 2.5; %162.64
    if (x < s) %hot [liquid]
        res = k1;
    else %cold [solid]
        res = k2;
    end;
end

