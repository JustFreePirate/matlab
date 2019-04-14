function res = K1( s, t )
    if (abs(s-t)<0.001)
        res = 1 / (2 + 3*t + t^2);
    else
        res = log((t+2)*(s+1)/((t+1)*(s+2))) / (s-t);
    end
end

