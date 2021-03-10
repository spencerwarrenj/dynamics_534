% This file illustrates the solution of a 2nd-order equation of motion
% using Euler's method. The method is applied to the solution of the
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
x0 = 1;			% Initial position (m)
v0 = 0;			% Initial velocity (m/s)

% Generate time vector
t = 0:dt:tf;
N = length(t);

% Initial conditions
y1(1) = x0;
y2(1) = v0;

% Step through each time step
for i = 1:N-1,
    
    % Calculate state derivatives
    y1dot = y2(i);
    y2dot = 1/m/(y1(i)^2+1)*(m*(omega^2-g)*y1(i)-m*y1(i)*y2(i)^2-k/2*y1(i)^3);
%     ydot = bead_wire_function(t(i), [y1(i) y2(i)]');
%     y1dot = ydot(1);
%     y2dot = ydot(2);
    
    % Calculate next value of position using Euler's method
	y1(i+1) = y1(i) + y1dot*dt;
	
	% Calculate next value of velocity using Euler's method
	y2(i+1) = y2(i) + y2dot*dt;

end

% Extract solution from y
x = y1;
v = y2;

% Plot solution
subplot(2,1,1);
plot(t,x,'b');
ylabel('x (m)');
subplot(2,1,2);
plot(t,v,'b');
xlabel('t (s)');
ylabel('v (m/s)');
