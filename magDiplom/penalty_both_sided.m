
N = 100; %x div
T = 1;
M = 100000; %time div    floor(T*N*N*2)
dx = 1 / N;
dt = T / M;
eps = 1e-4; % eps for penalty method
invEps = 1 / eps;

%if (x < s) hot [liquid]   u > 0
%if (x > s) %cold [solid]  u < 0
s0 = 0.5;
u0 = @(x) -1 * (x - 0.5);
I0 = @(x) -1 * (x*x/2 - 0.5*x); %% integrate of u0

kl = 0.56; %left k
kr = 2.5; %right k

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
    %solve for s    
    yls = vectorGetLerp(yl(k-1, :), s(k-1) / dx + 1); %for approx for -u_x(s(t)-, t)
    yrs = vectorGetLerp(yr(k-1, :), s(k-1) / dx + 1); %for approx for u_x(s(t), t)
    
    tk = t(k-1);

    dsdt = kl * invEps * yls + kr * invEps * yrs;
    %dsdt = 1 * invEps * yls + 1 * invEps * yrs;
    %dsdt = 0;
    
    s(k) = s(k-1) + dt * dsdt;
    
    
    %solve plain with straighforward Euler
    %inner points
    for i = 2:N
        yl_xx = (yl(k-1, i+1) - 2 * yl(k-1, i) + yl(k-1, i-1)) / dx^2;
        yr_xx = (yr(k-1, i+1) - 2 * yr(k-1, i) + yr(k-1, i-1)) / dx^2;
        
        %if (y(k-1, i) - ys > 0)
        if (x(i) <= s(k-1))
            xl = 0;
        else
            xl = 1;
        end
        
        if (x(i) >= s(k-1))
            xr = 0;
        else
            xr = 1;
        end

        yl(k, i) = yl(k-1, i) + dt * (kl * yl_xx - invEps * xl *(yl(k-1, i) - yls) - invEps * yls);
        yr(k, i) = yr(k-1, i) + dt * (kr * yr_xx - invEps * xr *(yr(k-1, i) - yrs) - invEps * yrs);
        
        %yl(k, i) = yl(k-1, i) + dt * (kl * yl_xx - invEps * min(yl(k-1, i) - yls, 0) - invEps * yls);
        %yr(k, i) = yr(k-1, i) + dt * (kr * yr_xx - invEps * max(yr(k-1, i) - yrs, 0) - invEps * yrs);
    end;
    %bounds
    yl(k, N+1) = 0;
    yr(k, 1) = 0;
    
    %NOTE: y_xx forward approximation
    yl(k, 1) = 2*yl(k, 2) - yl(k, 3);
    yr(k, N+1) = 2*yr(k, N) - yr(k, N-1);
    
    if (s(k) < 0 || s(k) > 1) 
        fprintf('s out of bound, s = %f, iter = %d \n', s(k), k);
        break;
    end;
end;


plot(t, s)

%disp(['s = ', s]);
fprintf('s = %f\n', s(M+1));

%for left: s last = 0.6151
%for right: s last = 0.2711  = 1 - 0.7289

p = 500
plot(yl(p, :))
hold on
plot(-fliplr(yr(p, :)))

