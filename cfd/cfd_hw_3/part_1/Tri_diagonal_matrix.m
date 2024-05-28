function x = Tri_diagonal_matrix(d, l, u, b)
    n = length(b);

    % Forward elimination
    for i = 2:n
        factor = l(i) / d(i - 1);
        d(i) = d(i) - factor * u(i - 1);
        b(i) = b(i) - factor * b(i - 1);
    end

    % Back substitution
    x(n) = b(n) / d(n);
    for i = n - 1:-1:1
        x(i) = (b(i) - u(i) * x(i + 1)) / d(i);
    end
end