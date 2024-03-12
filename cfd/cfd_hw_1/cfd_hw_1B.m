% Date  : 10.3.24
% Author: @eshnr7 (https://github.com/eshnr7)

clc; clear; close all;
%% cited:
%https://www.mathworks.com/matlabcentral/fileexchange/102189-solving-blasius-equation-using-rk-4-numerical-method

%% parameters 

% Defining ODEs
f_prime = @(f, g, h) g;
g_prime = @(f, g, h) h;
h_prime = @(f, g, h) -0.5*f*h;

% Defining initial conditions
f0 = 0;
g0 = 0;
h0 = 0.35; % Initial guess for h(0), it is obtained by trial and error

% defining η and Δη value, also stepsize
eta_0 = 0;
eta_max = 5;
d_eta = 0.001;
eta = eta_0:d_eta:eta_max;

F = zeros(length(eta), 3);
F(1, :) = [f0, g0, h0];
%% loop
% RK4 method
for i = 1:length(eta)-1
    eta_i = eta(i);
    f_i = F(i, 1);
    g_i = F(i, 2);
    h_i = F(i, 3);
    
    % RK4 steps
    k1 = d_eta * [f_prime(f_i, g_i, h_i), g_prime(f_i, g_i, h_i), h_prime(f_i, g_i, h_i)];
    k2 = d_eta * [f_prime(f_i + 0.5*k1(1), g_i + 0.5*k1(2), h_i + 0.5*k1(3)), ...
                  g_prime(f_i + 0.5*k1(1), g_i + 0.5*k1(2), h_i + 0.5*k1(3)), ...
                  h_prime(f_i + 0.5*k1(1), g_i + 0.5*k1(2), h_i + 0.5*k1(3))];
    k3 = d_eta * [f_prime(f_i + 0.5*k2(1), g_i + 0.5*k2(2), h_i + 0.5*k2(3)), ...
                  g_prime(f_i + 0.5*k2(1), g_i + 0.5*k2(2), h_i + 0.5*k2(3)), ...
                  h_prime(f_i + 0.5*k2(1), g_i + 0.5*k2(2), h_i + 0.5*k2(3))];
    k4 = d_eta * [f_prime(f_i + k3(1), g_i + k3(2), h_i + k3(3)), ...
                  g_prime(f_i + k3(1), g_i + k3(2), h_i + k3(3)), ...
                  h_prime(f_i + k3(1), g_i + k3(2), h_i + k3(3))];
              
    
    F(i+1, :) = F(i, :) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end
f = F(:, 1);
g = F(:, 2);
%% shear stress and skin friction coefficient
g_derivative = diff(transpose(g)) ./ diff(eta);
% g_derivative(1) = 0.35, this value is used numerically

x_1 = 0:.05:2;                  % defining x label for the shear_stress and c_f

shear_stress = zeros(length(x_1),1);
for i = 1:length(x_1)
    
    shear_stress(i) = 0.018*sqrt(1/x_1(i));
end

c_f = shear_stress/(.5*1.225*25);


%% Plotting
subplot(2,2,1)
plot(eta, f, 'LineWidth', 1); 
hold on;  
plot(eta, g, 'LineWidth', 1);             
title('Comparative plots of f(η) and g(η)');
xlabel('\eta');
legend('f(\eta)', 'g(\eta)'); 
hold off

subplot(2,2,2)
plot(g, eta)
xlabel("u/U∞")
ylabel("η")
legend("u(y)")
title("Comparative plots of velocity profile u(y) at x=1m")

subplot(2,2,3)
eta = eta(1:length(g));
plot(x_1, shear_stress)
xlabel("x")
ylabel("τ")
legend("τ(x)")
title("Shear Stress")

subplot(2,2,4)

plot(x_1, c_f)
xlabel("x")
ylabel("Cd_f")
legend("Cd_f(x)")
title("Skin friction coefficient")
