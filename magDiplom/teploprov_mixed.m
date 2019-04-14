% solve 
%   u_t = u_xx
%   u_x(0, t) = f(t)
%   u(1, t) = 0
%   u(x, 0) = u0(x)

%converges with much less iterations, N = 10, M can be 40 with w = 0.3,T=1
%u0 = @(x) -1 * (x - 1);
%f = @(t) 0;

N = 50; %x div
T = 1;
M = 1000; %time div   
dx = 1 / N;
dt = T / M;
w = 0.3; %relaxation factor
MAX_GZ_ITER = 200000;
GZ_EPS = dt * 1e-3;

u0 = @(x) (x - 1)^3/6;
f = @(t) t + 0.5;
exactU = @(x, t) (x-1) .* (t + (x-1).^2 ./ 6);

%u(k, :) - k-th layer [u( , t_k)]
u = zeros(M+1, N+1);
x = zeros(N+1, 1);

%fill first layers
% 1 layer of u and x
for i = 0:N
    xi = dx * i;
    x(i+1) = xi;
    u(1, i+1) = u0(xi);
end;


for k = 2:(M+1)
    tk = dt * (k - 1);
    
    %solve with over-relaxation method
    %wich is just Gauss-Zeidel but with interpolation between old solution
    %and new calculated solution with a factor of w
  
    for i = 1:N+1
        u(k, i) = u(k-1, i);
    end;
    for j = 1:MAX_GZ_ITER
        maxDiff = 0;
        %inner points
        for i = 2:N
            u_xx = (u(k, i+1) - 2 * u(k, i) + u(k, i-1)) / dx^2;
            prev = u(k, i);
            u(k, i) = (1 - w) * u(k, i) + w * (u(k-1, i) + dt * u_xx );
            maxDiff = max(maxDiff, abs(prev - u(k,i)));
        end;
        %bounds
        u(k, N+1) = 0;
        prev = u(k, 1);
        u(k, 1) = (1 - w) * u(k, 1) + w * (u(k, 2) - dx * f(tk));
        maxDiff = max(maxDiff, abs(prev - u(k, 1)));
        
        if (maxDiff < GZ_EPS)
            break;
        end;
    end;
    
    %inner points
    for i = 2:N
        u_xx = (u(k-1, i+1) - 2 * u(k-1, i) + u(k-1, i-1)) / dx^2;
        u(k, i) = (u(k-1, i) + dt * u_xx + u(k, i)) / 2;
    end;
    %bounds
    u(k, N+1) = 0;
    u(k, 1) = (u(k, 2) - dx * f(tk) + u(k, 1)) / 2;
    
end;

plot(x, u(M+1, :));
hold on;
plot(x, exactU(x, T))

