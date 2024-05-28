clear all, close all, clc

%% Information
% Date  : 24.3.24
% Author: @eshnr7 (https://github.com/eshnr7)

%% Parameters
rho     = 1000;             % Density
mu      = .001;             % Dynamic viscosity
nu      = mu / rho;         % Kinematic viscosity
delta_P = 5.12;             % Pressure drop
R       = 2.5 * 10^(-2);    % Radius
L       = 10;               % Length
Re      = 2000;             % Reynolds Number
T       = 650;              % Fully developed area, it has been found empirically

%% Definitions
nodes_number    = 50;                   % 5, 25, 50
dr              = 2*R/(nodes_number-1); % Defining dr
r               = -R:dr:R;              % r value distribution
dt              = 0.01;                 % Defining dt
time_number     = T / dt;               % Defining how many parts time will be divided into

U = zeros(time_number, nodes_number);   % Defining velocity

%% Boundary conditions
U(:, end)   = 0; % u(r = R, t) = 0
U(1, :)     = 0; % u(r, t = 0) = 0

%% Loop
for j = 1:time_number-1

    U(j, 1) = U(j+1, 1) - (dt/rho) * (delta_P / L) + ...
        (nu * dt / dr^2) * (U(j+1, 2) - 2 * U(j+1, 1) + U(j+1, end)) + ...
        (1 / (2 * r(1))) * (nu * dt / dr) * (U(j+1, 2) - U(j+1, end));

    U(j, end) = U(j+1, end) - (dt/rho) * (delta_P / L) + ...
        (nu * dt / dr^2) * (U(j+1, 1) - 2 * U(j+1, end) + U(j+1, end-1)) + ...
        (1 / (2 * r(end))) * (nu * dt / dr) * (U(j+1, 1) - U(j+1, end-1));
    
    for i = 2:nodes_number-1

        U(j+1, i) = U(j, i) - (dt/rho) * (delta_P / L) + ...
            (nu * dt / dr^2) * (U(j, i+1) - 2 * U(j, i) + U(j, i-1)) + ...
            (1 / (2 * r(i))) * (nu * dt / dr) * (U(j, i+1) - U(j, i-1));
    end
end
%% Wall Shear Stress Calculation

du_dr = (U(end, end) - U(end, end-1)) / dr; % calculating the du/dr for r=R
tau_w = mu * du_dr;
%% Flow Rate Calculation

radiues_matrix  = -R:dr:0; 
space_int       = radiues_matrix(1:nodes_number/2); 
Q               = 2*pi*trapz(space_int, U(end, 1:nodes_number/2).*space_int); 
%% Plotting

plot(-U(end, :), r, 'Color', 'Red', 'LineWidth', 1.25)
ylabel('Radius (m)');
xlabel('Velocity (m/s)');
title('Velocity Profile U(r)');

hold on;
plot([-min(U(end, :)) -max(U(end, :))], [0 0], 'black--');
hold off;
