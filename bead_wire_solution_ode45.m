% This file illustrates the solution of a 2nd-order equation of motion
% using ode45. The method is applied to the solution of the
% equation of motion of a bead on a rotating parabolic wire.
%
% Author:   Mark Colton
% Date:     2/10/21

clear;
clc;

% Make the model parameters global so that they can be accessed by the
% function in which the equations of motion are defined
global k m omega g;

% Define the model parameters and solution parameters
k = 50;			% Spring stiffness (N/m)
m = 0.5;		% Slider mass (kg)
omega = 0.5;	% Rotation rate (rad/s)
g = 9.81;		% Gravity acceleration (m/s^2)
dt = 0.001;		% Step size (s)
tf = 5;			% End time (s)
y0 = 1;			% Initial position (m)
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
[t,y] = ode45('bead_wire_function',t, y0);

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
