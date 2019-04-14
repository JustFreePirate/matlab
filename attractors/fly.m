syms V teta wz varteta x H Iz  real;
syms m S bA Psi g real;
syms P(V, H);
syms cx(varteta, teta, wz) cy(varteta, teta, wz) mz(varteta, teta, wz) q(V,H);
alpha = varteta - teta;




f1(V, teta, wz, varteta, x, H) = -g * sin(teta) + P(V, H)/m*cos(alpha) + q(V,H) * S/m *(-cx(varteta, teta, wz) * cos(alpha) - cy(varteta, teta, wz) * sin(alpha));
f2(V, teta, wz, varteta, x, H) = -g / V * cos(teta) + P(V, H)/(m * V) * sin(alpha) + q(V,H)*S/(m*V)*(-cx(varteta, teta, wz)*sin(alpha) +cy(varteta, teta, wz) * cos(alpha));
f3(V, teta, wz, varteta, x, H) = q(V,H) * S *bA * mz(varteta, teta, wz) / Iz;
f4(V, teta, wz, varteta, x, H) = wz;
f5(V, teta, wz, varteta, x, H) = V * cos(teta) * cos(Psi);
f6(V, teta, wz, varteta, x, H) = V * sin(teta);

F = [f1; f2; f3; f4; f5; f6];

j = jacobian(F);

disp(j);

