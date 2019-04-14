m = zeros(8);
m(1, 1) = 1;

for i = 1:8
    m(i, 1) = 1;
    m(1, i) = 1;
end;

for i = 2:8
    for j = 2:8
        m(i, j) = m(i-1,j) + m(i, j-1);
    end;
end;

disp(m);
