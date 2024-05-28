clear all, close all, clc

%% Information
% Date  : 29.4.24
% Author: @eshnr7 (https://github.com/eshnr7)

%% Parameters
dx = 0.05;
dt = 0.01;

K  = dt/dx; % coefficient

t_steps = [0 2 4 6];
Lx = 4;

Nx = Lx/dx + 1;
Nt = t_steps(end)/dt + 1;

%% Defining initial velocity distribution
u = zeros(Nx, Nt);
for i = 1:Nx
    x = (i-1)*dx;
    if x < 0.25
        u(i, 1) = 1;
    elseif x >= 0.25 && x <= 1.25
        u(i, 1) = 1.25 - x;
    else
        u(i, 1) = 0;
    end
end

%% Boundary conditions
u(1, :) = 1;
u(Nx, :) = 0;

%% Loop
for j = 1:Nt-1
    % Predictor step
    upred = u(:, j);
    for i = 2:Nx-1
        upred(i) = u(i, j) - K * 0.5 * (u(i+1, j)^2 - u(i, j)^2);
    end

    % Corrector step
    for i = 2:Nx-1
        u(i, j+1) = 0.5 * (u(i, j) + upred(i) - K * 0.5 * (upred(i)^2 - upred(i-1)^2));
    end
end

%% Plotting
figure;
for i = 1:length(t_steps)
    t_ind = t_steps(i)/dt + 1;
    plot(dx*(0:Nx-1), u(:, t_ind));
    hold on;
end
hold off;

xlabel('x (m)');
ylabel('u(x,t)');
title('Solution of Inviscid Burgers Equation using MacCormack Method');
legend('t=0', 't=2', 't=4', 't=6');
