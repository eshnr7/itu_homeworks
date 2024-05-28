clear all, close all, clc

% Please make sure that the "Tri_diagonal_matrix.m" file exists in the
% current directory.

%% Information
% Date  : 15.4.24
% Author: @eshnr7 (https://github.com/eshnr7) 

% Source Code: https://www.youtube.com/watch?v=NB0J3ENdbKQ&t=1277s
%% Parameters
time = 1000;
x    = 30;     
y    = 10;

% 30x10, 90x30, 180x60
n_y = 10;
n_x = 30;

d_t  = .01;         
d_x  = x/n_x;      
d_y  = y/n_y;      
%% Surface Temperatures
surface_1      = 18;      
surface_2      = 18;     
surface_3      = 18;     
surface_4      = 18;

surface_1_new  = 150;

max_temperature                = zeros(n_y, n_x);
max_temperature(end, :)        = surface_1_new;
max_temperature(:, end)        = surface_2;
max_temperature(1, :)          = surface_3;
max_temperature(:, 1)          = surface_4;

%% Defining Mesh
node_x  = (x / d_x) - 2;     
node_y  = (y / d_y) - 2;      
n       = node_x * node_y;               
node_t       = time / d_t;                
d       = d_t / (d_x^2);   % Diffusion number
%% Defining Tri-Diagonal Matrix

diagonal_A   =  2 * (1 + d) * ones(n, 1); % main diagonal
upper_dia_A  = -d * ones(n-1, 1);         % off diagonal
lower_dia_A  = -d * ones(n-1, 1);         % off diagonal

% Coefficient matrix of off diagonal
c     = zeros(3, 1);
c(1)  = d;
c(2)  = 2 * (1 - d);
c(3)  = d;

% Boundary Condition-1
b_c_1      = zeros(node_y, 1);
b_c_1(1)   = surface_3;
b_c_1(end) = surface_1_new;

% Boundary Condition-2
b_c_2      = zeros(node_x, 1);
b_c_2(1)   = surface_4;
b_c_2(end) = surface_2;

% Temperature
temp_distr_matrix  = zeros(n_y, n_x, node_t);
t   = 0;
%% loop
for count = 1:node_t
    for i = 1:node_x
        x_temp_distribution    = [max_temperature(1:end-2, i+1), max_temperature(2:end-1, i+1), max_temperature(3:end, i+1)];
        x_flux                 = x_temp_distribution * c + d * b_c_2(i);
        x_temp_update          = Tri_diagonal_matrix(diagonal_A, lower_dia_A, upper_dia_A, x_flux);
        
        max_temperature(2:end-1, i+1) = x_temp_update;
    end      
    for j = 1:node_y
        y_temp_distribution     = [max_temperature(j+1, 1:end-2); max_temperature(j+1, 2:end-1); max_temperature(j+1, 3:end)];
        c_resized               = repmat(c, 1, size(y_temp_distribution, 2));            
        y_flux                  = sum(y_temp_distribution .* c_resized, 1)' + d * b_c_1(j);
        y_temp_update           = Tri_diagonal_matrix(diagonal_A, lower_dia_A, upper_dia_A, y_flux);
        
        max_temperature(j+1, 2:end-1)  = y_temp_update';
    end
  
    temp_distr_matrix(:, :, count) = max_temperature;
    
    t = t + d_t;
end
%% Plotting

x_values = linspace(0, x, n_x); 
y_values = linspace(0, y, n_y); 

imagesc(x_values, y_values, temp_distr_matrix(:, :, end));
colorbar;
title('Attained Equilibrium State')
xlabel("X")
ylabel("Y")

%% The time required for the temperature of the center of the block to reach 30 degree

time_required = 0;

center_x_index = floor(n_x / 2) + 1;
center_y_index = floor(n_y / 2) + 1;

while temp_distr_matrix(center_y_index, center_x_index, time_required + 1) < 30
    time_required = time_required + 1;
    if time_required >= node_t
        disp("The temperature did not reach 30 degrees.");
        break;
    end
end

if time_required < node_t
    disp(['The temperature reaches 30 degrees at ', num2str(time_required * d_t), ' seconds']);
end
