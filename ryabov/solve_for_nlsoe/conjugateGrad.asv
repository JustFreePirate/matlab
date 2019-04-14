function conjugateGrad( eps )
    maxIter = 100000;
    iter = 0;
    x = [1 1 1]';
    f = F(x);
    err = norm(f);
    dx = 1;
    s = 0;

    while (err > eps && iter < maxIter)
        dx_prev = dx;
        s_prev = s;

        %compute J(x) (f should be computed)
        j = J(x);

        %compute grad and dx
        %grad F(v)^2/2 = J(v)'*F(v)
        grad = j' * f;
        dx = -grad;

        %compute beta
        %Fletcher-Reeves
        beta = dot(dx, dx) / dot(dx_prev, dx_prev);

        %Polak-Ribiere
        %beta = dot(dx, dx - dx_prev) / dot(dx_prev, dx_prev);

        %compute direction s
        s = dx + beta * s_prev;

        %compute alpha (line search for argmin{alpha} F^2(x + alpha * s))
        v = j * s;
        alpha = -dot(f, v) / dot(v, v);

        %update position
        x = x + alpha * s;

        f = F(x);
        err = norm(f);
        iter = iter + 1;
    end;

    disp('x = ');
    disp(x);

    disp('F(x) = ');
    disp(F(x));

    disp('iter = ');
    disp(iter);
end