function res = foo( t, x )
    M = 1;
    m = 2;
    l = 1;
    g = 9.8;
    tau = [9.8*3 0]';
    
    %x = [q1 q2 q1' q2']^t = [u a u' a']^t
    u = x(1);
    a = x(2);
    du = x(3);
    da = x(4);
    Mq = [  m+M         m*l*cos(a)
            m*l*cos(a)  m*(1+(l*sin(a))^2)];
    
    v = [-m*l*sin(a)*da^2 + m*g + M*g
         m*l^2*sin(2*a)*da^2-m*l*sin(a)*da*du-m*(l*l*da*da*sin(2*a)/2-l*du*da*sin(a))+m*g*l*cos(a)];
    r34 = Mq^(-1)*(tau - v);
    res = [ du
            da
            r34(1)
            r34(2)];
end

