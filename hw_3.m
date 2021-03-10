%hw_3 %newtonian dynamics

%% Problem 3.1
clear;
t = 3;
theta = pi/20*cos(2*t);
theta_dot = -pi/10*sin(2*t);
w1 = .5;
w3 = 7;
R = mry(theta)*mrx(-50*pi/180);
omega_disk_rel_o_in_Fprime = [w3;theta_dot;0];
omega_disk_rel_inert_in_bigF = R'*(omega_disk_rel_o_in_Fprime)+[0;0;w1]

%% Problem 3.2
clear;
t = pi; %sec
y = 40;  %cm
y_dot = -30; %cm/s
y_ddot = -4; %cm/s^2
w = .2; %rad/s
theta = pi/6*sin(2*t);
theta_dot = pi/3*cos(2*t);
theta_ddot = -2*pi/3*sin(2*t);
R = mrx(theta);
omega_frame = R*[0;0;w]+[theta_dot;0;0];
r = [0;y;0];
r_dot = [0;y_dot;0]+cross(omega_frame,r)
r_dot_ana = [-y*w*cos(theta);y_dot;theta_dot*y]
r_ddot = [-y_dot*w*cos(theta)+y*w*sin(theta)*theta_dot; ...
    y_ddot; ...
    theta_ddot*y+theta_dot*y_dot]+cross(omega_frame,r_dot)
syms theta theta_dot theta_ddot w t y y_dot y_ddot
R = mrx(theta);
omega_frame = R*[0;0;w]+[theta_dot;0;0];
r = [0;y;0];
r_dot = [0;y_dot;0]+cross(omega_frame,r)
r_ddot = [-y_dot*w*cos(theta)+y*w*sin(theta)*theta_dot; ...
    y_ddot; ...
    theta_ddot*y+theta_dot*y_dot]+cross(omega_frame,r_dot)
cross(omega_frame,r_dot)

%% Problem 3.4
clear;
clc;

% Make the model parameters global so that they can be accessed by the
% function in which the equations of motion are defined
global mu w g;

% Define the model parameters and solution parameters
mu = .5; %coefficient of friction
w = 10; %omega (rad/s)
g = 9.81;		% Gravity acceleration (m/s^2)
dt = 0.001;		% Step size (s)
tf = .8;			% End time (s)
y0 = .05;		% Initial position (m)
v0 = 0;			% Initial velocity (m/s)

% Generate time vector
t = 0:dt:tf;
N = length(t);

% Initial conditions
y0 = [y0 v0];

% Call 'ode45' to solve ODE.  'ode45' calls the function
% 'bead_wire_function' repeatedly, which returns the two 
% derivatives at each time step.  'ode45' uses the returned
% derivatives to calculate the solution at each time step.
% The solution is returned in the variable y, which contains
% both y1 (x) and y2 (v).
[t,y] = ode45(@hw3_4_function,t, y0);

% Extract solution from y
x = y(:,1);		% Position
v = y(:,2);		% Velocity

% Plot solution
subplot(2,1,1);
plot(t,x,'r');
ylabel('x (m)');
subplot(2,1,2);
plot(t,v,'r');
xlabel('t (s)');
ylabel('v (m/s)');

%% Problem 3.5
clear; clc; close;
global m k R;
m = .75; %kg
k = 25; %N/m
R = 1; %m
g = 9.81;		% Gravity acceleration (m/s^2)
dt = 0.0001;		% Step size (s)
tf = acos((2*pi-5)/(-5));			% End time (s)
y0 = 0;		% Initial position (m)
v0 = 0;			% Initial velocity (m/s)

% Generate time vector
t = 0:dt:tf;
N = length(t);

% Initial conditions
y0 = [y0 v0];

% Call 'ode45' to solve ODE.  'ode45' calls the function
% 'bead_wire_function' repeatedly, which returns the two 
% derivatives at each time step.  'ode45' uses the returned
% derivatives to calculate the solution at each time step.
% The solution is returned in the variable y, which contains
% both y1 (x) and y2 (v).
[t,y] = ode45(@hw3_5_function,t, y0);


theta = -5.*cos(t)+5;
w = 5*sin(t);
w_dot = 5*cos(t);

% Extract solution from y
s = y(:,1);		% Position
s_dot = y(:,2);		% Velocity

% N = m.*sqrt((2.*s_dot.*w+s.*w_dot-R.*w.^2).^2+g.^2);
N = m.*(2.*s_dot.*w+s.*w_dot-R.*w.^2);

figure(1)
subplot(2,1,1)
plot(t,s)
xlabel('time (s)')
ylabel('s (m)')
subplot(2,1,2)
plot(t,N)
xlabel('time (s)')
ylabel('Normal Force (N)')

figure(2)
subplot(2,1,1)
plot(theta,s)
xlabel('\theta (rad)')
ylabel('s (m)')
subplot(2,1,2)
plot(theta,N)
xlabel('\theta (rad)')
ylabel('Normal Force (N)')

max(abs(s))
max(abs(N))

%% Problem 3.6
clear; close; clc;

syms theta
solve(1==(2-theta)*cos(theta)+sin(theta),theta)

%% functions
function R = mrz(theta)
    R = [cos(theta), sin(theta), 0;
        -sin(theta), cos(theta), 0;
        0, 0, 1];
end
function R = mrx(theta)
    R = [1,0,0;
        0,cos(theta), sin(theta);
        0,-sin(theta),cos(theta)];
end
function R = mry(theta)
    R = [cos(theta),0,-sin(theta);
        0,1, 0;
        sin(theta),0,cos(theta)];
end