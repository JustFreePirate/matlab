syms alpha w teta
syms a1 a2 a3 a4

A = [-a1 1 0; -a2 -a3 0; 0 1 0];
C = eye(3);
D = [0 0 0];
B = [0 -a4 0];

sys = ss(A, C, B, D);

