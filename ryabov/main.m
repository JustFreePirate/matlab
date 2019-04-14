%z(t) = t(t-1)
%z(0) = z(1) = 1

K = @(s,t) 1 ./ (s+t+1);
%K = @(s, t) 1 ./ (2 - s*t);
u = @(s) -s-1.5+(s.^2+3.*s+2).*log((s+2)./(s+1));
%u = @(s) (s.^2 - 4.*s.*(1 + log(2)) + log(256) + 4.*(-2 + s).*log(2 - s))./(2.*s.^3);
p = @(s) 1 + s*100;
r = @(s) 1;
alpha = 1e-10;

n = 25;
h = 1 / n;
U = zeros(n+1, 1);
U(1) = 0;
U(n+1) = 0;
S = zeros(1, n+1);
M = zeros(n+1, n+1);
for k = 1:n+1
    S(k) = (k-1)*h;
end;

M(1,1) = 1;
M(n+1, n+1) = 1;


for k = 2:n 
    sk = S(k);
    U(k) = integral(@(t) K(t, sk).*u(t), 0, 1);
    pl = p(sk - h/2);
    pr = p(sk + h/2);
    rk = r(sk);
    M(k, k-1) = alpha * (-pl) / h^2;
    M(k, k) = alpha * ((pl + pr) / h^2 + rk);
    M(k, k+1) = alpha * (-pr) / h^2;
    
    M(k, 1) = M(k, 1) + h / 2 * K1(sk, 0);
    M(k, n+1) = M(k, n+1) + h / 2 * K1(sk, 1);
    for j = 2:n
        tj = S(j);
        M(k, j) = M(k, j) + h * K1(sk, tj); 
    end;
end

Z = M^(-1)*U;
zReal = zeros(n+1, 1);
N = zeros(n+1, 1);
for k = 1:n+1
    sk = S(k);
    zReal(k) = sk*(sk-1);
    for j = 2:n
        tj = S(j);
        N(k) = N(k) + h * K(sk, tj) * Z(j);
    end
    N(k) = abs(N(k) - u(sk));
end;

plot(S, Z, S, zReal, S, N);
legend('zk','z real', 'int-u');


