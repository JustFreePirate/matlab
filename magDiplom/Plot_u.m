function Plot_u(x, u, s, us, k )
    dx = x(2) - x(1);
    lastIndex = size(u, 2);
    leftIndex = floor(s(k) / dx) + 1;
    rightIndex = ceil(s(k) / dx) + 1;
    if (leftIndex == rightIndex) 
        leftIndex = leftIndex - 1;
        rightIndex = rightIndex + 1;
        %заменяем серединку в u
    end;
    
    xm = [x(1:leftIndex); s(k); x(rightIndex:lastIndex)];
    um = [u(k, 1:leftIndex), us(k), u(k, rightIndex:lastIndex)];
    plot(xm, um);
    hold on;
    plot(s(k), us(k), 'r*'); 
end

