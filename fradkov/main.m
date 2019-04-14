[T, Y] = ode45(@(t,y) foo(t,y), [0 5], [0; 0; 0; 0]);
plot(T, Y(:,1), '-', T, Y(:,2), '-');
title('Modeling');
xlabel('Time t');
ylabel('Solution y');
legend('u','alpha');