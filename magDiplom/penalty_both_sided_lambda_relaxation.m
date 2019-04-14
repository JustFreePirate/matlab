
N = 50; %x div
T = 1;
M = 3000; %time div    floor(T*N*N*2)
dx = 1 / N;
dt = T / M;
eps = 1e-4; % eps for penalty method
invEps = 1 / eps;
lambda = 0; %-0.006

w = 0.2; %relaxation factor
MAX_GZ_ITER = 200000;
GZ_EPS = dt * 1e-3;

%if (x < s) hot [liquid]   u > 0
%if (x > s) %cold [solid]  u < 0
s0 = 0.5;
u0 = @(x) -2 * (x - 0.5);
I0 = @(x) -2 * (x*x/2 - 0.5*x); %% integrate of u0

kl = 1; %left k  0.56
kr = 2; %right k  2.5

%kl = 1;
%kr = 1;

%ul - u in (0, s(t))
%ur - u in (s(t), 1)
%yl(x, t) = int from x to 1 ul(ksi, t) dKsi
%yr(x, t) = int from 0 to x ur(ksi, t) dKsi

yl = zeros(M+1, N+1);
yr = zeros(M+1, N+1);
x = zeros(N+1, 1);
s = zeros(M+1, 1);
t = zeros(M+1, 1);
us = zeros(M+1, 1);


%fill first layers
% 1 layer of u and x
for i = 0:N
    xi = dx * i;
    x(i+1) = xi;
    if (xi < s0)
        yl0 = I0(s0) - I0(xi);
        yr0 = 0;
    else
        yl0 = 0;
        yr0 = I0(xi) - I0(s0);
    end
    
    yl(1, i+1) = yl0;
    yr(1, i+1) = yr0;
end;
s(1) = s0;

for k = 1:M+1
    t(k) = dt * (k - 1);
end;


for k = 2:(M+1)
    s(k) = s(k-1);
    for i = 1:N+1
        yl(k, i) = yl(k-1, i);
        yr(k, i) = yr(k-1, i);
    end;
    
    for j = 1:MAX_GZ_ITER
        maxDiff = 0;
        
        %solve for s    
        yls = vectorGetLerp(yl(k, :), s(k) / dx + 1); %for approx for -u_x(s(t)-, t)
        yrs = vectorGetLerp(yr(k, :), s(k) / dx + 1); %for approx for u_x(s(t), t)

        dsdt = kl * invEps * yls + kr * invEps * yrs;
        us(k) = -lambda * dsdt;
        
        prev_s = s(k);
        new_s = s(k-1) + dt * dsdt;
        s(k) = (1 - w) * prev_s + w * new_s;
        
        if (s(k) < 0 || s(k) > 1) 
            break;
        end;
        
        maxDiff = max(maxDiff, abs(prev_s - new_s));
        
        if (k == 2)
            ddsddt = (s(k) - 2 * s(k-1) + s(k-1)) / (dt*dt);
        else
            ddsddt = (s(k) - 2 * s(k-1) + s(k-2)) / (dt*dt);
        end;

        %solve plain with straighforward Euler
        %inner points
        for i = 2:N
            yl_xx = (yl(k, i+1) - 2 * yl(k, i) + yl(k, i-1)) / dx^2;
            yr_xx = (yr(k, i+1) - 2 * yr(k, i) + yr(k, i-1)) / dx^2;

            if (x(i) <= s(k))
                Al = invEps * yls - lambda * ddsddt * (s(k) - x(i));
            else
                Al = invEps * yl(k, i);
            end

            if (x(i) >= s(k))
                Ar = invEps * yrs - lambda * ddsddt * (x(i) - s(k));
            else
                Ar = invEps * yr(k, i);
            end

            prev_yl = yl(k, i);
            new_yl = yl(k-1, i) + dt * (kl * yl_xx - Al);
            yl(k, i) = (1 - w) * prev_yl + w * new_yl;
            
            prev_yr = yr(k, i);
            new_yr = yr(k-1, i) + dt * (kr * yr_xx - Ar);
            yr(k, i) = (1 - w) * prev_yr + w * new_yr;
            
            maxDiff = max(maxDiff, abs(prev_yl - new_yl));
            maxDiff = max(maxDiff, abs(prev_yr - new_yr));
        end;
        %bounds
        yl(k, N+1) = 0;
        yr(k, 1) = 0;

        %NOTE: y_xx forward approximation
        prev_yl = yl(k, 1);
        new_yl = 2*yl(k, 2) - yl(k, 3);
        yl(k, 1) = (1 - w) * prev_yl + w * new_yl;
        
        prev_yr = yr(k, N+1);
        new_yr = 2*yr(k, N) - yr(k, N-1);
        yr(k, N+1) = (1 - w) * prev_yr + w * new_yr;
        
        maxDiff = max(maxDiff, abs(prev_yl - new_yl));
        maxDiff = max(maxDiff, abs(prev_yr - new_yr));
        
        if (maxDiff < GZ_EPS)
            break;
        end;
    end
    
    if (s(k) < 0 || s(k) > 1) 
        fprintf('s out of bound, s = %f, iter = %d \n', s(k), k);
        break;
    end;
end;

fprintf('lambda = %f\n', lambda);
plot(t, s)

%plot(t, us); %display temp in s


%disp(['s = ', s]);
fprintf('s = %f\n', s(M+1));

%for left: s last = 0.6151
%for right: s last = 0.2711  = 1 - 0.7289

% p = 500
% plot(yl(p, :))
% hold on
% plot(-fliplr(yr(p, :)))

