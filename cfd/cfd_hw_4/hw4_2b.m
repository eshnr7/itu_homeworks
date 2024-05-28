clear all, close all, clc
%% Information
% Date  : 29.4.24
% Author: @eshnr7 (https://github.com/eshnr7)

%% Defining the function u(x, t)
u = @(x, t) sin(2*pi*x) .* (cos(2*pi*t) + sin(2*pi*t));

x = 0:0.01:1;
t_values = [0.1, 0.2, 0.3];

%% Plot
hold on;
for i = 1:length(t_values)
    plot(x, u(x, t_values(i)), 'DisplayName', ['time = ' num2str(t_values(i))]);
end

xlabel('x');
ylabel('u(x, t)');
title('Analitic Solution');
legend('show');
hold off;