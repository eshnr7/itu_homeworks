function x = jacobi_method(a1, b1, tol, max_iter, x0)
    x = x0;
    for iter = 1:max_iter
        x_new = zeros(size(b1));
        for i = 1:length(b1)
            
            x_new(i) = (b1(i) - a1(i,:) * x + a1(i,i) * x(i)) / a1(i,i);
        end
        if norm(x_new - x) < tol
            break;
        end
        x = x_new;
    end
end
