clear all, close all  , clc
%________________________________________
% Parameters:

m = 14000;              %kg
lamda = 1.3;
g = 9.80665;            %m/s^2
v_i = 36*(5/18);        %m/s
v_flat = 60*(5/18);     %m/s
theta_i = 7;            %degree
theta_o = 0;            %degree
A = 0.8;                %m^2
c = 500;                %J/KgK
h = 60;                 %W/m^2K
T_a = 20+273;           %Kelvin
m_b = 30;               %kg
f_r = 0.01;
A_p = 5;                %m^2
C_D = 0.5;
row = 1.2;              %kg/m^3
a_target = 6;           %m/s^2
a_i = 0; %speed is constant for first scenario
%________________________________________
% 1.A:

F_i = lamda*m*a_i;
F_G = m*g*sind(theta_i);
F_R = f_r*m*g*cosd(theta_i);
% -F_B-F_R+F_G+F_i = 0 => F_B = -F_R+F_G+F_i

F_B = -F_R+F_G+F_i; % that's equal to 15369 N
%________________________________________
% 1.B:

F_B_front = F_B*0.6;
F_B_rear = F_B*0.4;

P_B_front = (F_B_front)*v_i; %that's equal to 92215 
P_B_rear = (F_B_rear)*v_i;   %that's equal to 61476
P_B_total = P_B_front + P_B_rear; %that's equal to 153690

%________________________________________
% 1.C:
t = 0:10:600; %defining time interval
%front disc:
T_front = 1921.14 + 293 - (1921.14./(exp(3.2e-3 * t))); %function found in report

subplot(2, 2, 1)
plot(t, T_front)

%rear disc:
T_rear = 1280.75 + 293 - (1280.75./(exp(3.2e-3 * t)));%function found in report

hold on
subplot(2, 2, 1)
plot(t, T_rear) 
hold off

title("Brake disc temperature-Analytically;1")
xlabel("Time (second)")
ylabel("Temperature")
legend("Front Disc", "Rear Disc")

%________________________________________
% 1.D:

h_step_size = 10;
t_final = 6000/10;
t(1) = 0;
T(1) = 293;

%----------------------------------------
% for front disc
f = @(t, T) (P_B_front+ h*A*(T_a-T))/(m_b*c);

for i = 1: ceil(t_final/h_step_size)
    
    t(i+1) = t(i) + h_step_size;

    k_1 = f(t(i), T(i));
    k_2 = f(t(i) + 0.5*h_step_size, T(i) + 0.5*k_1*h_step_size);
    k_3 = f(t(i) + 0.5*h_step_size, T(i) + 0.5*k_2*h_step_size);
    k_4 = f(t(i) + h_step_size, T(i) + k_3*h_step_size);

    T(i+1) = T(i) + (h_step_size/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
end

subplot(2, 2, 2)
plot(t,T)
%----------------------------------------
% for rear disc

f = @(t, T) (P_B_rear+ h*A*(T_a-T))/(m_b*c);

for i = 1: ceil(t_final/h_step_size)
    
    t(i+1) = t(i) + h_step_size;

    k_1 = f(t(i), T(i));
    k_2 = f(t(i) + 0.5*h_step_size, T(i) + 0.5*k_1*h_step_size);
    k_3 = f(t(i) + 0.5*h_step_size, T(i) + 0.5*k_2*h_step_size);
    k_4 = f(t(i) + h_step_size, T(i) + k_3*h_step_size);

    T(i+1) = T(i) + (h_step_size/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
end
hold on
subplot(2, 2, 2)
plot(t,T)

title("Brake disc temperature-Runge-Kutta;1")
xlabel("Time (second)")
ylabel("Temperature")
legend("Front Disc", "Rear Disc")
hold off
%________________________________________
% 2.A:
% forces are writen according to the new values

F_G_2 = m*g*sind(theta_o);
F_R_2 = f_r*m*g*cosd(theta_o);

F_i_2 = lamda*m*a_target;

t_final_2 = v_flat/a_target; %new t_final
h_step_size_2 = 0.01; % new step_size


f = @(t, T, P_B_front_2) (P_B_front_2 + h * A * (T_a - T)) / (m_b * c);

for i = 1:ceil(t_final_2 / h_step_size_2)

    t(i+1) = t(i) + h_step_size_2;

    v_2 = v_flat - a_target * t(i);

    F_AD = 0.5 * row * C_D * A_p * v_2^2;
    F_B_2 = -F_R_2 + F_G_2 + F_i_2 - F_AD;

    F_B_front_2 = F_B_2 * 0.6;
    P_B_front_2 = F_B_front_2 * v_2;

    k_1 = f(t(i), T(i), P_B_front_2);
    k_2 = f(t(i) + 0.5 * h_step_size_2, T(i) + 0.5 * k_1 * h_step_size_2, P_B_front_2);
    k_3 = f(t(i) + 0.5 * h_step_size_2, T(i) + 0.5 * k_2 * h_step_size_2, P_B_front_2);
    k_4 = f(t(i) + h_step_size_2, T(i) + k_3 * h_step_size_2, P_B_front_2);

    T(i+1) = T(i) + (h_step_size_2 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
end

subplot(2, 2, 3)
plot(t,T)

%----------------------------------------
% for rear disc

f = @(t, T, P_B_rear_2) (P_B_rear_2 + h * A * (T_a - T)) / (m_b * c);

for i = 1:ceil(t_final_2 / h_step_size_2)
    
    t(i+1) = t(i) + h_step_size_2;

    v_2 = v_flat - a_target * t(i);

    F_AD = 0.5 * row * C_D * A_p * v_2^2;
    F_B_2 = -F_R_2 + F_G_2 + F_i_2 - F_AD;

    F_B_rear_2 = F_B_2*0.4;
    P_B_rear_2 = (F_B_rear_2)*v_2;

    k_1 = f(t(i), T(i), P_B_rear_2);
    k_2 = f(t(i) + 0.5 * h_step_size_2, T(i) + 0.5 * k_1 * h_step_size_2, P_B_rear_2);
    k_3 = f(t(i) + 0.5 * h_step_size_2, T(i) + 0.5 * k_2 * h_step_size_2, P_B_rear_2);
    k_4 = f(t(i) + h_step_size_2, T(i) + k_3 * h_step_size_2, P_B_rear_2);

    T(i+1) = T(i) + (h_step_size_2 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
end
hold on
subplot(2, 2, 3)
plot(t,T)
hold off

title("Brake disc temperature-Runge-Kutta;2")
xlabel("Time (second)")
ylabel("Temperature")
legend("Front Disc", "Rear Disc")
