function x = SOR_method_Func(A, b, omega, tol, max_iter)

    x = zeros(size(b));
    n = length(b);
    iter = 0;
    residual = inf;

    while residual > tol && iter < max_iter
        x_old = x;
        for i = 1:n
            sigma = A(i,1:i-1) * x(1:i-1) + A(i,i+1:n) * x_old(i+1:n);
            x(i) = (1 - omega) * x_old(i) + (omega / A(i,i)) * (b(i) - sigma);
        end 

        residual = norm(A * x - b);
        iter = iter + 1;
    end
end