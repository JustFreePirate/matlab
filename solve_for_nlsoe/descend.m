function descend( eps )
    maxIter = 100000;
    iter = 0;
    x = [1 1 1]';
    f = F(x);
    err = norm(f);
    while (err > eps && iter < maxIter)
        j = J(x);
        grad = j' * f;
        v = -j * grad;
        alpha = -dot(f, v) / dot(v, v);
        x = x - alpha * grad;
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