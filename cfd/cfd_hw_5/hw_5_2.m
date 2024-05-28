% Date  : 28.5.24
% Author: @eshnr7 (https://github.com/eshnr7)

clear all, close all, clc
%% Parameters
mu      = 1.0;      %When Re=1, then mu=rho
rho     = 1.0;
Re      = 1.0; 
L       = 1.0;
h       = 0.02;
dt      = 0.0001;   % Timestep size, selected by trial and error to ensure stability condition
t_end   = 15;
beta    = 1.6;      % Overrelaxation factor for SOR method

%% Main code in stability condition control
d = 2 * (mu / rho) * dt / h^2;

if d > 0.5
    disp('The stability criterion d > 0.5 inequality was not met.')
else
    
    nx = round(L / h) + 1;
    ny = round(L / h) + 1;
    
    u = zeros(ny, nx); 
    v = zeros(ny, nx); 
    psi = zeros(ny, nx); 
    omega = zeros(ny, nx); 
    
    %% Boundary conditions
    u(:, 1) = 0; 
    v(:, 1) = 0; 
    u(:, nx) = 0; 
    v(:, nx) = 0; 
    u(ny, :) = 0; 
    v(ny, :) = 0; 
    u(1, :) = 1; 
    v(1, :) = 0; 
    
%% Loop
    psi_new = psi; 
    omega_new = omega; 
    iter_count = 0;

    for t = 0:dt:t_end
        iter_count = iter_count + 1;
        
        %% Updating stream function using SOR method
        for i = 2:nx-1
            for j = 2:ny-1
                psi_new(j, i) = 0.25 * beta * ...
                    (psi(j, i+1) + psi_new(j, i-1) + ...
                    psi(j+1, i) + psi_new(j-1, i) + ...
                    h^2 * omega_new(j, i)) + ...
                    (1 - beta) * psi(j, i);
            end
        end
        
        %% Updating vorticity boundaries
        for i = 1:nx
            omega(1, i) = 2 * (psi(1, i) - psi(2, i)) / h^2 + 2 * u(1, i) / h;
            omega(ny, i) = 2 * (psi(ny, i) - psi(ny-1, i)) / h^2 - 2 * u(ny, i) / h;
        end
        for j = 1:ny
            omega(j, 1) = 2 * (psi(j, 1) - psi(j, 2)) / h^2 - 2 * v(j, 1) / h;
            omega(j, nx) = 2 * (psi(j, nx) - psi(j, nx-1)) / h^2 + 2 * v(j, nx) / h;
        end
        
        %% Updating vorticity using FTCS
        for i = 2:nx-1
            for j = 2:ny-1
                omega_new(j, i) = omega(j, i) + dt * ...
                    ((mu / rho) * ...
                    ((omega(j, i+1) - 2 * omega(j, i) + omega(j, i-1)) / h^2 + ...
                    (omega(j+1, i) - 2 * omega(j, i) + omega(j-1, i)) / h^2) - ...
                    (u(j, i) * (omega(j, i+1) - omega(j, i-1)) / (2 * h) - ...
                    v(j, i) * (omega(j-1, i) - omega(j+1, i)) / (2 * h)));
            end
        end
        
        psi = psi_new;
        omega = omega_new;
    end
    
    %% Calculating u and v from stream function
    for i = 2:nx-1
        for j = 2:ny-1
            u(j, i) = -(psi(j-1, i) - psi(j+1, i)) / (2 * h);
            v(j, i) = (psi(j, i+1) - psi(j, i-1)) / (2 * h);
        end
    end
    
% Define grid coordinates
[x, y] = meshgrid(0:h:L, 0:h:L);

%% Plotting
    
    % Stream function
    figure;
    contourf(psi_new, 'ShowText', 'on', 'LineWidth', 2);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Stream Function');
    set(gca, 'YDir', 'reverse');
    
    % Vorticity
    figure;
    contourf(omega_new, 'ShowText', 'on', 'LineWidth', 2);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Vorticity');
    set(gca, 'YDir', 'reverse');
    
    % Velocity vectors 
    figure;
    quiver(u, v);
    xlim([0 nx]);
    ylim([0 ny]);
    xlabel('x');
    ylabel('y');
    title('Velocity Vectors');
    set(gca, 'YDir', 'reverse');
    
    % u velocity
    figure;
    contourf(u, 'ShowText', 'on', 'LineWidth', 2);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('u Velocity');
    set(gca, 'YDir', 'reverse');
    
    % v velocity 
    figure;
    contourf(v, 'ShowText', 'on', 'LineWidth', 2);
    colorbar;
    xlabel('x');
    ylabel('y');
    title('v Velocity');
    set(gca, 'YDir', 'reverse');



end
