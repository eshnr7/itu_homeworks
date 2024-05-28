% Date  : 27.5.24
% Author: @eshnr7 (https://github.com/eshnr7)
% benefited from the first pseudocode (p. 223, lecture slides)
clear all, close all, clc
%% Parameters 

N1 = 101;
M1 = 101;
L = 1.0;
H = L / (N1 - 1); % also equal to K
Re_list = [1, 5, 10, 50]; % List of Reynolds 

%% Grid setup
x = linspace(0, L, N1);
y = linspace(0, L, M1);
[X, Y] = meshgrid(x, y);

%% Defining U and V velocities
U = zeros(N1, M1); 
V = zeros(N1, M1);

%% Boundary conditions
U(1, :) = 1;        % Top boundary moving right
U(end, :) = 0;      % Bottom boundary stationary
U(:, 1) = 0;        % Left boundary stationary
U(:, end) = 0;      % Right boundary stationary
V(1, :) = 0;        % Top boundary
V(end, :) = 0;      % Bottom boundary
V(:, 1) = 0;        % Left boundary
V(:, end) = 0;      % Right boundary


figure;

%% Loop 
for idx = 1:length(Re_list)
    RE = Re_list(idx);
    A1 = 0.5 * RE * H;
    G = 1; % considered as H = K so H^2/K^2 = 1 
    A2 = 2.0 * (1 + G);
    GR = sqrt(G);

    for iter = 1:1000
        for I = 2:N1-1
            for J = 2:M1-1
                % Update U
                B1_u = 1.0 - A1 * U(I, J);
                B3_u = 2.0 - B1_u;
                B2_u = G - A1 * GR^2 * V(I, J);
                B4_u = 2.0 * G - B2_u;

                U(I, J) = (U(I+1, J) * B1_u + U(I-1, J) * B3_u + U(I, J+1) * B2_u + U(I, J-1) * B4_u) / A2;

                % Update V
                B1_v = 1.0 - A1 * V(I, J);
                B3_v = 2.0 - B1_v;
                B2_v = G - A1 * GR^2 * U(I, J);
                B4_v = 2.0 * G - B2_v;

                V(I, J) = (V(I+1, J) * B1_v + V(I-1, J) * B3_v + V(I, J+1) * B2_v + V(I, J-1) * B4_v) / A2;
            end
        end
%% Checking for convergence
        if mod(iter, 100) == 0
            fprintf('Iteration %d for Re=%d\n', iter, RE);
        end
    end
%% Plotting
    subplot(2, 2, idx);
    contourf(X, Y, U, 20);
    colorbar;
    hold on;
    %quiver(X, Y, U, V, 'k');
    title(sprintf('Velocity Field for Re = %d', RE));
    xlabel('x');
    ylabel('y');
    %axis equal tight;
    hold off;
end