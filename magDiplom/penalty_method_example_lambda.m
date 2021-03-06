%solve for
%u_t = u_xx
%u(0, t) = t + 0.5
%u_x(s(t),t)= -3 * s_t + 2t/sqrt(3-2t)
%u(x,0) = (x^2/2 - 2x + 0.5)^+
%exact solution: 
%   u(x,t) = (x^2/2 - 2x + 0.5 + t)^+
%   s(t) = 2 - sqrt(3-2t)

%y(x,t) = int from x to R u(ksi, t)dKsi

N = 50; %x div
T = 1;
M = 10000; %time div    floor(T*N*N*2)
dx = 1 / N;
dt = T / M;
eps = 1e-4; % eps for penalty method
invEps = 1 / eps;
lambda = 1e-5;

f = @(t) t + 0.5;
exactU = @(x, t) max(0, x.^2 ./ 2 - 2 .* x + 0.5 + t);

h = @(t) (2*t/sqrt(3-2*t)) / (3);
s0 = 2 - sqrt(3);
%I = @(x, t) x.^3./6 - x.^2 + x./2 + x.*t;
%exactY = @(x, t) max(0, I(exactS(t), t) - I(x, t));

%u(k, :) - k-th layer [u( , t_k)]
y = zeros(M+1, N+1);
x = zeros(N+1, 1);
s = zeros(M+1, 1);
t = zeros(M+1, 1);

dsdt_log = zeros(M+1, 1);

%fill first layers
% 1 layer of u and x
for i = 0:N
    xi = dx * i;
    x(i+1) = xi;
    y(1, i+1) = y0(xi);
end;
s(1) = s0;

for k = 1:M+1
    t(k) = dt * (k - 1);
end;



for k = 2:3%(M+1)
    %solve for s
    ys = vectorGetLerp(y(k-1, :), s(k-1) / dx + 1);
    
    tk = t(k-1);
    %ys = exactY(s(k-1), tk);
    
    %ys = sqrt(3-2*tk);
    %ys = 0;
    dsdt = (1/3 * invEps * ys + h(tk)) / (1 - invEps * lambda * (1 - s(k-1)));
    
    %dsdt = sqrt(3-2*tk)/3 +h(tk);
    
    dsdt_log(k) = ys;
    s(k) = s(k-1) + dt * dsdt;
    
    %solve plain with straighforward Euler
    %inner points
    for i = 2:N
        y_xx = (y(k-1, i+1) - 2 * y(k-1, i) + y(k-1, i-1)) / dx^2;
        
        if (x(i) <= s(k-1))
            Aeps = invEps * ys + invEps * lambda * dsdt * (1 - s(k-1));
        else
            Aeps = invEps * y(k-1, i) + invEps * lambda * dsdt * (1 - x(i));
        end
        
        y(k, i) = y(k-1, i) + dt * (y_xx - Aeps);
    end;
    %bounds
    y(k, N+1) = y(k, N) + lambda * dsdt * dx ;
    y(k, 1) = y(k, 2) + dx * f(t(k));
    %y(k, 1) = exactY(0, t(k));
    
    if (s(k) <= 0 || s(k) >= 1) 
        fprintf('s out of bound, s = %f, iter = %d \n', s(k), k);
        break;
    end;
end;

% plot(t, s)
% hold on
% plot(t, exactS(t))

% yy = zeros(N+1, 1);
% for i = 1:N+1
%     yy(i) = exactY(x(i), 0.5);
% end
% 
% plot(x, y(M / 2,:))
% hold on
% plot(x, yy)



function res = exactS(t) 
    res = 2 - sqrt(3 - 2.*t);
end

function res = exactY(x, t) 
    I = @(x, t) x.^3./6 - x.^2 + x./2 + x.*t;
    s = exactS(t);
    if x < s
        res = I(s, t) - I(x, t);
    else
        res = 0;
    end;
end

function res = y0(x)
    r = 2 - sqrt(3);
    if (x < r)
        int = @(x) x^3/6 - x^2 + x/2;
        res = int(r) - int(x);
    else
        res = 0;
    end;
end