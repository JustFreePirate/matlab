% solve 
%   u_t = u_xx
%   u_x(0, t) = f(t)
%   u(1, t) = 0
%   u(x, 0) = u0(x)

%converges with N = 10, T = 1, M = 200 (T * N * N * 2),
%u0 = @(x) -1 * (x - 1);
%f = @(t) 0;
%but diverges with M = 100;

N = 50; %x div
T = 1;
M = 5000; %time div    floor(T*N*N*2)
dx = 1 / N;
dt = T / M;

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
    %solve plain with straighforward Euler
    %inner points
    for i = 2:N
        u_xx = (u(k-1, i+1) - 2 * u(k-1, i) + u(k-1, i-1)) / dx^2;
        u(k, i) = u(k-1, i) + dt * u_xx;
    end;
    %bounds
    u(k, N+1) = 0;
    u(k, 1) = u(k, 2) - dx * f(tk);
end;

plot(x, u(M+1, :));
hold on;
plot(x, exactU(x, T))


