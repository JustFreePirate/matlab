function plot_uu(x, yl, yr, s, us, k )
    [M, N] = size(yl);
    dx = x(2) - x(1);
    lastIndex = N;
    leftIndex = floor(s(k) / dx) + 1;
    rightIndex = ceil(s(k) / dx) + 1;
    if (leftIndex == rightIndex) 
        leftIndex = leftIndex - 1;
        rightIndex = rightIndex + 1;
        %�������� ��������� � u
    end;
    
    %convert to temp yl
    ul = zeros(1, N);
    ur = zeros(1, N);
    
    for i = 1:N-1
        ul(i) = -(yl(k, i+1) - yl(k, i)) / dx  + us(k);
    end;
    for i = N:-1:2
        ur(i) = (yr(k, i) - yr(k, i-1)) / dx  + us(k);
    end;
    

    %todo: add us(k)
%approximation with h^2
%     ul(1) = -(yl(k, 2) - yl(k, 1)) / dx;
%     for i = 2:N-1
%         ul(i) = -(yl(k, i+1) - yl(k, i-1)) / (2*dx);
%     end;
%     ur(N) = (yr(k, N) - yr(k, N-1)) / dx;
%     for i = N-1:-1:2
%         ur(i) = (yr(k, i+1) - yr(k, i-1)) / (2*dx);
%     end;
    
    
    xm = [x(1:leftIndex); s(k); x(rightIndex:lastIndex)];
    um = [ul(1:leftIndex), us(k), ur(rightIndex:lastIndex)];
    plot(xm, um);
    hold on;
    plot(s(k), us(k), 'r*'); 
end

