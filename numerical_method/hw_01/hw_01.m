clear all, close all, clc
%-------------------------------------------------
% parameters:

m_round = 0.28*10^(-3); %kg
d = 0.006;              %m
T = 16 + 273.15;        %K
fi = 58/100;
P = 1024*100;           %pa
X_0 = 0;
Y_0 = 1.5;
dt = 0.01;
Ar = (pi*d^2)/4;        %m^2
Cd = 0.545;             % determined by trial and error
g = 9.81;               %m/s^2
%-------------------------------------------------
% Calculating air density:

%thermodynamic table-a4 for T = 16, P_saturated = 1832.4 Pa
e =  fi*1832.4; 
% to calculate air_density 
air_density = (P-e)/(287*T); %here 287 is R, table-a1

%-------------------------------------------------
% a) Curve fitting: 

x0 = 1; x1 = 6; x2 = 11; x3 = 16;
f_x0 = 108.83; f_x1 = 88.84; f_x2 = 76.44; f_x3 = 64.26667;

xs = [x0, x1, x2, x3];
f_x = [f_x0, f_x1, f_x2, f_x3];

p_0 = polyfit(xs, f_x, 2);

x_values = linspace(0,16,100);
y_values = polyval(p_0,x_values);


equation = p_0(1)+"x^2\t" + p_0(2)+"x\t" + p_0(3); %that's the equation used experimental data
fprintf(equation)

%-------------------------------------------------
% c) Euler Forward Formula

%initialize variables
v_x = 112.6343; %initial value
v_y = 0;
t = 0;
x = X_0;
y = Y_0;
times = t;
velocities_x = v_x;
velocities_y = v_y;
positions_x = x;
positions_y = y;

%kinematic equations
while y > 0
    
    t = t + dt;

    v_y = v_y +(((0.5*air_density*Cd*Ar*v_y^2)-(m_round*g))/m_round)*dt;
    v_x = v_x + (-0.5*air_density*Cd*Ar*v_x^2/m_round)*dt;
    
    y = y + v_y*dt;
    x = x + v_x*dt;

    times = [times t];
    
    velocities_x = [velocities_x v_x];
    positions_x = [positions_x x];

    velocities_y = [velocities_y v_y];
    positions_y = [positions_y y];
end

%-------------------------------------------------
%Graphics:

subplot(2,1,1)
plot(xs,f_x,'*',x_values,y_values)
xlabel('X(m)')
ylabel('V(x) m/s')

hold on
plot(positions_x, velocities_x)
title("CD:   "+Cd)
legend('Experimantal data','Curve fitting', 'Modeled graphic')

subplot(2,1,2)
plot(positions_x, positions_y)
xlabel('X(m)')
ylabel('Y(m)')
legend("Modeled graphic")

