% Date  : 27.5.24
% Author: @eshnr7 (https://github.com/eshnr7)
% benefited from 
% the first pseudocode  (p. 223, lecture slides)
% the second pseudocode (p. 231, lecture slides)
clear all, close all, clc
%% Parameters
N1 = 101;
M1 = 101; 
L = 1.0; 
H = L / (N1 - 1); % also equal to K
Re_list = [100, 1000, 10000]; % List of Reynolds

%% Grid setup
x = linspace(0, L, N1);
y = linspace(0, L, M1);
[X, Y] = meshgrid(x, y);

%% Defining U and V velocities
U = zeros(N1, M1);
V = zeros(N1, M1); 

%% Boundary conditions
U(1, :) = 1;    % Top boundary moving right
U(end, :) = 0;  % Bottom boundary stationary
U(:, 1) = 0;    % Left boundary stationary
U(:, end) = 0;  % Right boundary stationary
V(1, :) = 0;    % Top boundary
V(end, :) = 0;  % Bottom boundary
V(:, 1) = 0;    % Left boundary
V(:, end) = 0;  % Right boundary

%% Loop 
for RE = Re_list
    A1 = 0.5 * RE * H;
    G = 1; % considered as H = K so H^2/K^2 = 1 
    A2 = 2.0 * (1 + G);
    GR = sqrt(G);

    F1 = H*RE;
    F2 = H*RE; % normally there is K instead of H, but they are same

    for iter = 1:1000
        for I = 2:N1-1
            for J = 2:M1-1
                %% Update U, 
                % here, the second pseudocode
                X1 = U(I,J);
                X2 = V(I,J);
                X3 = 1 + F1*abs(X1);
                X4 = G*(1.0 + F2*abs(X2));
             
                if X1 > 0.0
                    B1_u = 1.0;
                    B3_u = X3;
                else
                    B1_u = X3;
                    B3_u = 1.0;
                end

                if X2 > 0.0
                    B2_u = G;
                    B4_u = X4;
                else
                    B2_u = X4;
                    B4_u = G;
                end

                B0_u = X3+1.0+X4+G;

                U(I, J) = (U(I+1, J) * B1_u + U(I-1, J) * B3_u + U(I, J+1) * B2_u + U(I, J-1) * B4_u) / B0_u;

                % Update V
                X5 = 1 + F1*abs(X2);
                X6 = G*(1.0 + F2*abs(X1));

                if X2 > 0.0
                    B1_v = 1.0;
                    B3_v = X5;
                else
                    B1_v = X5;
                    B3_v = 1.0;
                end

                if X1 > 0.0
                    B2_v = G;
                    B4_v = X6;
                else
                    B2_v = X6;
                    B4_v = G;
                end

                B0_v = X5+1.0+X6+G;              

                V(I, J) = (V(I+1, J) * B1_v + V(I-1, J) * B3_v + V(I, J+1) * B2_v + V(I, J-1) * B4_v) / B0_v;
            end
        end
%% Checking for convergence
        if mod(iter, 100) == 0
            fprintf('Iteration %d for Re=%d\n', iter, RE);
        end
    end

%% Plotting
    figure;
    contourf(X, Y, U, 20);
    colorbar;
    hold on;
    quiver(X, Y, U, V, 'k');
    title(sprintf('Velocity Field for Re = %d', RE));
    xlabel('x');
    ylabel('y');
    axis equal tight;
    hold off;
end
