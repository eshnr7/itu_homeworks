clear all, close all, clc
% Please make sure that the "SOR_method_Func", "jacobi_method", and 
% "gauss_siedel_method" files exists in the current directory.
%% Information
% Date  : 15.4.24
% Author: @eshnr7 (https://github.com/eshnr7) 

% Source Code: https://engineering-stream.com/matlab_codes_links/Laplace_Equation_2D_Dirichlet_BCs.html
%% Parameters
x = 30;
y = 10;

% 30x10, 90x30, 180x60
n_x = 30;
n_y = 10;

d_x = x/n_x; 
d_y = y/n_y;

%% Surface Temperature
T1 = 18;
T2= 18;
T3 = 150;
T4 = 18;

%% Defining Mesh
m=((x/d_x) - 1); 
n= ((y/d_y) - 1); 
r=n*m; 

%% Creating a matrix
for i=2:r+1
    for j= 2:r+1
        if i == j
            a(i,j) = -4;
        elseif (i == j+1)&&((round((i-2)/n)) ~= (i-2)/n)
            a(i,j)=1;
        elseif i == j-1 &&((round((i-2)/n)) ~= (j-2)/n)
            a(i,j)=1;
        elseif i == j+n
            a(i,j) =1;   
        elseif i == j-n
            a(i,j)=1;
        else
            a(i,j) = 0;
        end
    end
end

a1=a(2:r+1,2:r+1);

%% Creating b vector
for i=2:m+1
    for j = 2:n+1
        if (i == 2) &&(j == 2)
            d(i,j) =- (T1+T4);
        elseif (j == 2) && ((i> 2) && (i<(m+1)))
            d(i,j) = -T1;
        elseif (j == 2) && (i == (m+1))
            d(i,j) =- (T1+T2);
        elseif (i == (m+1)) && ((j> 2) && (j< (n+1)))
            d(i,j)= -T2;
        elseif (i == (m+1)) && (j == (n+1))
            d(i,j) =- (T2+T3);
        elseif (j == (n+1)) && ((i> 2) && (i<(m+1)))
            d(i,j)= -T3;
        elseif (i == 2) && (j == (n+1))
            d(i,j) =- (T3+T4);
        elseif (i == 2) && ((j > 2) && (j < (n+1)))
            d(i,j) =- T4;
        else
            d(i,j) = 0;
        end
    end
end

matrix = d(2:m+1, 2:n+1);
k = 2;
for i = 2:m+1
    for j = 2:n+1
        b(i+j-k)=d(i,j);
    end
    k=k-(n-1);
end

b1=b(2:r+1)';

%% Prepering the Methods

omega    = 1.6;             % Relaxation parameter, it is given in Homework file.
tol      = 1e-6;            % Tolerance for convergence
max_iter = 1000;            % Maximum number of iterations
x0       = zeros(size(b1)); % using in the Jacobi Method

%% Solving the system using Jacobi method

%solution = jacobi_method(a1, b1, tol, max_iter, x0);

%% Solving the system using Gauss-Seidel method

%solution = gauss_siedel_method(a1, b1, tol, max_iter, x0);

%% Solving the system using Gauss-Seidel SOR method

solution = SOR_method_Func(a1, b1, omega, tol, max_iter);

%% Creating T matrix
k=2;
for i = 2:m+1
    for j = 2:n+1
        T(i,j)=solution(i+j-k-1);
    end
    k=k-(n-1);
end

%% Plotting
T(2:m+1,2:n+1);
T12= T(2:m+1, 2:n+1)';

x_values = linspace(0, x, n_x); 
y_values = linspace(0, y, n_y); 

imagesc(x_values,y_values,T12)
colorbar
title('Attained Equilibrium State')
xlabel("X")
ylabel("Y")

%% Checking for convergence
delta_T = inf;
iter = 0;
while delta_T > tol && iter < max_iter

    old_solution = solution;
    solution = gauss_siedel_method(a1, b1, tol, 1, solution);
    
    delta_T = abs(solution(ceil(m/2) * n + ceil(n/2)) - old_solution(ceil(m/2) * n + ceil(n/2)));
    
    iter = iter + 1;
end

fprintf('Temperature at center node converged to: %.6f\n', solution(ceil(m/2) * n + ceil(n/2)))
fprintf('Number of iterations: %d\n', iter)