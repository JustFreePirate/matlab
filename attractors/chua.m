clear;
alpha = 8.4562;
beta = 12.0732;
gamma = 0.0052;
m0_num = -0.1768;
m1_num = -1.1468;

syms a b g p w m1 m0 t k x;
assume([a b g p w m1 m0 t k x], 'real');
P = [-a*(1+m1) a 0; 1 -1 1; 0 -b -g];
q = [-a; 0; 0];
r = [1 0 0];
W(p) = r * inv(P - p*eye(3)) * q;
roots = solve(imag(W(t*1i))==0, t);
roots = roots(2:3);%% only positive roots
f(a,b,g) = roots;
w_nums = double(f(alpha, beta, gamma));

f(a,b,g,m1) = - 1./real(W(1i*w_nums));
k_nums = double(f(alpha, beta, gamma, m1_num)); 

%w = [2.0392 3.2453]
%k = [0.2099 0.9597]

% P0 = P + k*q*r;
% W0(p) = r*inv(P0 - p*eye(3))*q;
% A = [0 -w 0; w 0 0; 0 0 -d];
% b = [b1; b2; 1];
% c = [1 0 -h];
% params: b1 b2 b3 A3 c3 w
% Wa(p) = c*inv(A - p*eye(3))*b; 

d = (a+w^2-b+1+g+g^2)/(1+g);
h = a*(g+b-(1+g)*d+d^2)/(w^2+d^2);
b1 = a*(g+b-w^2-(1+g)*d)/(w^2+d^2);
f(a, b, g, w) = b1;
b1_nums = double(f(alpha, beta, gamma, w_nums));

s11 = 1; s12 = 0; s13 = -h;
s21 = m1 + 1 + k; s22 = -w/a;
s23 = -h*(a*(m1+1+k)-d)/a;
s31 = (a*(m1+k)-w^2)/a;
s32 = -(a*(b+g)*(m1+k)+a*b-g*w^2)/(a*w);
s33 = h*(a*(m1+k)*(d-1)+d*(1+a-d))/a;
S = [s11 s12 s13; s21 s22 s23; s31 s32 s33];

%a0 = 9.4287 for root1

syms x
assume(x, 'real');
F(x, t) = 0.5*(x+sin(x)*cos(x))+t*sin(x);
I(t) = F(acos(-t), t) - F(0, t) - (F(2*pi-acos(-t), t) - F(acos(-t), t)) + F(2*pi, t) - F(2*pi-acos(-t), t);

J(a) = 0.5*a*(I(1/a)-I(-1/a));
Fi(a) = ((m0-m1)*J(a)-a*k*pi)/w;
f(a, m0, m1, k, w) = Fi(a);

%solve: Fi(a)==0

a1 = solve(f(a, m0_num, m1_num, k_nums(1), w_nums(1))==0, a);
%a1 = 5.8561450862573606726561747674329


a2 = solve(f(a, m0_num, m1_num, k_nums(2), w_nums(2))==0, a);
%a2 = 1.0449210798936908291958771822724

%dFi(a) = diff(Fi, a)
%dFi(a1) = -0.3201 < 0
%dFi(a2) = -0.3319 < 0
%f(a, m0, m1, w, k) = dFi
%dFia1 = double(f(a1, m0_num, m1_num, w_nums(1), k_nums(1)))
%dFia2 = double(f(a2, m0_num, m1_num, w_nums(2), k_nums(2)))
%a2 не подходит, т.к. b1(w2) < 0

Y = [a1 0 0]';
X = S * Y;
f(a, m1, k, w) = X;
X0 = double(f(alpha, m1_num, k_nums(1), w_nums(1)));
X0 = -X0;
%X0 = 
% 5.8561450862573606726561747674329 
% 0.36933157824675905223886754751202
% -8.3665361683320104844752434709817

% solve   X' = f(t, X) = P0*X+eps*q*phi(r*x)
P0 = P + k*q*r;
f(a, b, g, k, m1) = P0;
P0 = double(f(alpha, beta, gamma, k_nums(1), m1_num));
% P0 =
% 
%    -0.5333    8.4562         0
%     1.0000   -1.0000    1.0000
%          0  -12.0732   -0.0052

f(a) = q;
q = double(f(alpha));
%phi(x) = (m0 - m1) * sat(x) - k*x
phi = @(x) (m0_num - m1_num) * 0.5 * (abs(x+1)-abs(x-1)) - k_nums(1)*x;

f = @(t, x, e) P0*x + e*q*phi(r*x);

T = 700;
eps = 0.1;
deltaEps = 0.1;
[~, x] = ode45(@(t, x) f(t, x, eps), [0 T], X0);

for n = 1:9
    X0 = x(length(x), :)';
    eps = eps + deltaEps;
    [~, x] = ode45(@(t, x) f(t, x, eps), [0 T], X0);
end;

plot3(x(:, 2), x(:, 1), x(:, 3),'-')
grid on




