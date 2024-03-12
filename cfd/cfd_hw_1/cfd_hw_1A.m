% Date  : 10.3.24
% Author: @eshnr7 (https://github.com/eshnr7)

clear all, close all, clc
%% parameters

% NOTE: 
% Dynamic viscosity, density, and U_infinity are not define as
% parameters, their numerical values are used direct.
%
% the matrix should be this format: [m]*[g_values] = [l]

N = 100;
m = zeros(N);                           % defining a N*N matrix
l = zeros(N, 1);                        % defining l matrix
h = 0:.05:5;                            % defining η and Δη value
r = @(x) x - (1/2)*sqrt(pi)*erf(x);     % this is integral of g funtion, it is taken from Wolfram Alpha
%% loops and calculations
for i = 1:N
    x = h(i);
    f = x - (1/2)*sqrt(pi)*erf(x);      % this equation is the ingtegral of "g = 1 - exp(-x.^2)"
    a = -4;                             % this is taken from using centred fdd numerical method
    if i == N-1
        b = 2 + (f/2)*(h(i+1) - h(i));
        l(N) = -b;                      % please see the left side of the matrix equation on page 5 of the report
    else
        b = 2 + (f/2)*(h(i+1) - h(i));  % this is taken from using centred fdd numerical method
        c = 2 - (f/2)*(h(i+1) - h(i));  % this is taken from using centred fdd numerical method
    end
    for j = 1:N
        if i == j+1
            m(i, j) = c;
        elseif i == j
            m(i, j) = a;
        elseif i == j-1
            m(i, j) = b;
        end
    end

end

g_function = m \ l;             % [g] = [m]^-1*[l]
h = h(1:length(g_function));
r_values = zeros(size(h));
for i = 1:length(h)
    r_values(i) = r(h(i));
end

g_derivative = diff(transpose(g_function)) ./ diff(h);
% g_derivative(1) = 0.3827, this value is used numerically

x_1 = 0:.05:2;                  % defining x label for the shear_stress and c_f

shear_stress = zeros(length(x_1),1);
for i = 1:length(x_1)
    
    shear_stress(i) = 0.02*sqrt(1/x_1(i));
end

c_f = shear_stress/(.5*1.225*25);
%% plotting

subplot(2,2,1)
plot(h, g_function)

hold on
plot(h, r_values)
hold off

legend("f'", "f")
xlabel('η')
ylabel('g and f values')
title('Comparative plots of f(η) and g(η)')

subplot(2,2,2)
plot(g_function, h)
xlabel("u/U∞")
ylabel("η")
legend("u(y)")
title("Comparative plots of velocity profile u(y) at x=1m")

subplot(2,2,3)
h = h(1:length(g_derivative));
plot(x_1, shear_stress)
xlabel("x")
ylabel("τ")
legend("τ(x)")
title("Shear Stress")

subplot(2,2,4)
h = h(1:length(g_derivative));
plot(x_1, c_f)
xlabel("x")
ylabel("Cd_f")
legend("Cd_f(x)")
title("Skin friction coefficient")