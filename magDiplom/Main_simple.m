%function stab = Main_simple( l )
clear;

l = 0.3;

sPoints = [];

%table of lambda stability
stabN = [10,     50,     100,    150,      200,        250,        300,         350];
stabL = [0.398,  0.403,  0.409,  0.415,    0.4168,     0.4158,     0.4149,      0.4144]; %*dx

N = 100; %x div
T = 2;
M = floor(T*N*N*6); %time div
dx = 1 / N;
dt = T / M;
% stable lambda is increasing func of N and decreasing of M
lambda = dx  * l; 

%if (x < s) hot [liquid]   u > 0
%if (x > s) %cold [solid]  u < 0
s0 = 0.5;
u0 = @(x) -1 * (x - 0.5);
w0 = @(x) 0;
w1 = @(x) 1*sin(pi * x);
k1 = 0.56; %28.21
k2 = 2.5; %162.64

%stability if dt / dx^2 <= 0.5
if (2*T*N*N*k2 > M) 
    disp('stability criteria is not met');
end;

%u(k, :) - k-th layer [u( , t_k)]
u = zeros(M+1, N+1);
w = zeros(M+1, N+1);
s = zeros(M+1, 1);
us = zeros(M+1, 1); %u(s(k), tk);
x = zeros(N+1, 1);

%fill first layers
% 1 layer of u
% 1 layer of w from w0 
for i = 0:N
    xi = dx * i;
    x(i+1) = xi;
    u(1, i+1) = u0(xi);
    w(1, i+1) = w0(xi);
end;
s(1) = s0;
us(1) = u0(s(1));

for k = 2:(M+1)
%for k = 2:3
    %calc w_k from w_k-1, w_k-2, u_k-1, s_k-1
    if (k == 2)
        %calc 2 layer of w from w1 and 1 layer
        %interior
        for i = 2:N
            xi = x(i);
            w(2, i) = w(1, i) + dt * w1(xi);
        end;
    else
        %interior
        for i = 2:N
            xi = x(i);
            sigma = Sigma(xi, u(k-1, i), s(k-1)); %sigma_k-1
            w_xx = (w(k-1, i+1) - 2 * w(k-1, i) + w(k-1, i-1)) / dx^2; %w_xx(x_i, t_k-1)
            den = 1/dt^2 + sigma / (2*dt);
            num = w_xx + sigma / (2*dt) * w(k-2, i) - (-2 * w(k-1, i) + w(k-2, i)) / dt^2;
            w(k, i) = num / den;
        end;
    end;
    
    %need to use u_s to calculate proper temp u
    %so we calculate where is s(k-1)
%     leftIndex = floor(s(k-1) / dx) + 1;
%     rightIndex = ceil(s(k-1) / dx) + 1;
%     if (leftIndex == rightIndex)
%         u(k-1, leftIndex) = us(k-1);
%         leftIndex = leftIndex - 1;
%         rightIndex = rightIndex + 1;
%     end;
    

    %calc s first (aka semi-implicit euler)
    sIndex = round(s(k-1) / dx) + 1;
    u_xl = (u(k-1, sIndex) - u(k-1, sIndex - 1)) / dx;
    u_xr = (u(k-1, sIndex + 1) - u(k-1, sIndex)) / dx;
    
    ds = -k1 * u_xl + k2 * u_xr;
    s(k) = s(k-1) + dt * ds;
    
    
    %calc layer u_k from u_k-1, s_k-1, w_k, w_k-1
    %interior
    for i = 2:N
        xi = x(i);
        sigma = Sigma(xi, u(k-1, i), s(k-1)); %sigma_k-1
        w_t = (w(k, i) - w(k-1, i)) / dt;
        %w_t = 0;
       
        %kappa = Kappa(xi, s(k-1));
        kappa = Kappa(xi, s(k));
        %kappa = Kappa(xi, round(s(k-1) / dx) * dx);
        %kappa = 2.5;
        
%         uc = u(k-1, i);
%         if (i == leftIndex) 
%             h = s(k-1) - x(i);
%             ur = us(k-1);
%             ul = lerp(u(k-1, i), u(k-1, i-1), h / dx);
%         elseif (i == rightIndex)
%             h = x(i) - s(k-1);
%             ur = lerp(u(k-1, i), u(k-1, i+1), h / dx);
%             ul = us(k-1);
%         else
%             h = dx;
%             ul = u(k-1, i-1);
%             ur = u(k-1, i+1);
%         end;
        
        %h = dx;
        %ul = u(k-1, i-1);
        %ur = u(k-1, i+1);
        
        u_xx = (u(k-1, i+1) - 2 * u(k-1, i) + u(k-1, i-1)) / dx^2;
        
        %u_xx = (ur - 2 * uc + ul) / h^2;
        
        u(k, i) = u(k-1, i) + dt * (kappa * u_xx + sigma * w_t^2);
    end;
    %boundary (u_x = 0 through boundary)
    u(k, 1) = u(k, 2);
    u(k, N+1) = u(k, N);
    
    %calc s_k from u_k-1, us_k-1, s_k-1
%     eps = 1e-5;
%     leftIndex = floor(s(k-1) / dx - eps) + 1;
%     rightIndex = ceil(s(k-1) / dx + eps) + 1;
%     u_xl = (us(k-1) - u(k-1, leftIndex)) / (s(k-1) - x(leftIndex));
%     u_xr = (u(k-1, rightIndex) - us(k-1)) / (x(rightIndex) - s(k-1));
    
   
    
    
    %us(k) = 0;
    
    us(k) = -lambda * ds * rand(1);
    
    if (s(k) <= 0 || s(k) >= 1)
        break;
    end
    
    
    %!!!!!round S!!!!!
    %s(k) = round(s(k)/dx)*dx;
    sIndexPrev = sIndex;
    sIndex = round(s(k) / dx) + 1;
    u(k, sIndex) = us(k);
    
    if (sIndexPrev ~= sIndex)
        sPoints(end+1) = k;
    end
    
    if (sIndex == 1 || sIndex == N+1) 
        break;
    end;
    
end;

if (k == (M+1)) 
    disp('soshlos');
else
    disp('ne soshols');
end;

% plot(x, u(3, :));
% hold on;
% plot(s(3), us(3), 'r*'); 
% 
% plot(x, u(40, :));
% plot(s(40), us(40), 'r*'); 

%end
