%solve for
%u_t = u_xx
%u(0, t) = t + 0.5
%u_x(s(t),t)= -3 * s_t + 2t/sqrt(3-2t)
%u(x,0) = (x^2/2 - 2x + 0.5)^+
%exact solution: 
%   u(x,t) = (x^2/2 - 2x + 0.5 + t)^+
%   s(t) = 2 - sqrt(3-2t)

%y(x,t) = int from x to R u(ksi, t)dKsi

N = 200; %x div
T = 1;
M = 300000; %time div    floor(T*N*N*2)
dx = 1 / N;
dt = T / M;
eps = 1e-5; % eps for penalty method
invEps = 1 / eps;

f = @(t) t + 0.5;
exactU = @(x, t) max(0, x.^2 ./ 2 - 2 .* x + 0.5 + t);
u0 = @(x) exactU(x, 0);

h = @(t) (2*t/sqrt(3-2*t)) / (3);
s0 = 2 - sqrt(3);
%I = @(x, t) x.^3./6 - x.^2 + x./2 + x.*t;
%exactY = @(x, t) max(0, I(exactS(t), t) - I(x, t));

%u(k, :) - k-th layer [u( , t_k)]
u = zeros(M+1, N+1);
x = zeros(N+1, 1);
s = zeros(M+1, 1);
t = zeros(M+1, 1);

dsdt_log = zeros(M+1, 1);

%fill first layers
% 1 layer of u and x
for i = 0:N
    xi = dx * i;
    x(i+1) = xi;
    u(1, i+1) = u0(xi);
end;
s(1) = s0;

for k = 1:M+1
    t(k) = dt * (k - 1);
end;



for k = 2:(M+1)
    %solve for s
    %calc int from s to 1 u(x)
    sIdx = s(k-1) / dx + 1;
    us = vectorGetLerp(u(k-1, :), sIdx);
    r = ceil(sIdx);
    %int_u = dx * (r - sIdx) * (us + u(k-1, r)) / 2;
    int_u = dx * (0 + u(k-1, r)) / 2;
    
    %integral with trapezium
    %summ inner points
    for i = r+1:N
        int_u = int_u + u(k-1, i) * dx ;
    end
    %summ boundaries
    if (r < N+1)
        int_u = int_u + (u(k-1, r) + u(k-1, N+1)) * dx / 2 ;
    end
    
    tk = t(k-1);
    dsdt = 1/3 * invEps * int_u + h(tk);

    s(k) = s(k-1) + dt * dsdt;
    %s(k) = exactS(t(k));
    
    %solve plain with straighforward Euler
    %inner points
    for i = 2:N
        u_xx = (u(k-1, i+1) - 2*u(k-1, i) + u(k-1, i-1)) / dx^2;
        
        %if (y(k-1, i) - ys > 0)
        if (x(i) <= s(k-1))
            xe = 0;
        else
            xe = 1;
        end

        u(k, i) = u(k-1, i) + dt * (u_xx - invEps * xe * u(k-1, i));
    end;
    %bounds
    u(k, N+1) = 0;
    u(k, 1) = f(t(k));

    if (s(k) < 0 || s(k) > 1) 
        fprintf('s out of bound, s = %f, iter = %d \n', s(k), k);
        break;
    end;
end;

plot(t, s)
hold on
plot(t, exactS(t))

% tt = 0.2;
% plot(x, u(ceil(M * tt), :));
% hold on;
% plot(x, exactU(x, tt));

%plot(t, dsdt_log);


function res = exactS(t) 
    res = 2 - sqrt(3 - 2.*t);
end
