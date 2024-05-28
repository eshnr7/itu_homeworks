clear all, close all, clc
%% Information
% Date  : 29.4.24
% Author: @eshnr7 (https://github.com/eshnr7)

%% Parameters
L = 1; % 0<x<1 
T = 1; 

dx = 0.01; 
dt = 0.01;

Nx = L/dx + 1; 
Nt = T/dt;

u = zeros(Nx,Nt);

%% Applying IC
x = linspace(0, L, Nx);

int_du_dt = -cos(2*pi*x);
u(:,1)    = int_du_dt;
u(:,1)    = sin(2*pi*x); % u(x, t=0) = sin(2*pi*x)
 
%% Applying BC
u(1,:)   = 0; % u(x=0, t) = 0
u(end,:) = 0; % u(x=1, t) = 0

%% Loop
for j = 1:Nt-1
    for i = 2:Nx-1
        u(i,j+1) = u(i,j)-(dt/(2*dx))*(u(i+1,j)-u(i-1,j))+...
            (dt^2/(2*dx^2))*(u(i+1,j)-2*u(i,j)+u(i-1,j));
    end
end

%% Plot

% Defining time range
t_values = [0.1, 0.2, 0.3];

hold on;
for i = 1:length(t_values)
    plot(x, u(:, round(t_values(i)/dt) + 1), 'DisplayName', ['time = ' num2str(t_values(i))]);
end

xlabel('x');
ylabel('u(x, t)');
title('Lax-Wendroff method');
legend('show');
hold off;
