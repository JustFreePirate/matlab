function [ y ] = F( v )
    x = v(1);
    y = v(2);
    z = v(3);
    y1 = x^2+y^2+z^2-1;
    y2 = 2*x^2+y^2-z^2;
    y3 = x^2 - exp(y)+1;
    y = [y1 y2 y3]';
end

