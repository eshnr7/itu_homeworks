function x = gauss_siedel_method(a, b, tol, max_iter, x0)

    x = x0;
    
    for iter = 1:max_iter
        x_new = x; 
        for i = 1:length(b)
            sum_ = 0;
            for j = 1:length(b)
                if j ~= i
                    sum_ = sum_ + a(i,j) * x_new(j);
                end
            end

            x_new(i) = (b(i) - sum_) / a(i,i);
        end
 
        if norm(x_new - x) < tol
            break;
        end
        x = x_new; 
    end
end
