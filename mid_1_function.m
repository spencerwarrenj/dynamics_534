% mid_1_function
% This is the function that 'ode45' calls to get the function derivatives
% at each time step.  'ode45' passes the current time and states (x and v,
% which are contained in y), and the function returns the current state
% derivatives (xdot and vdot, which are contained in ydot).
%
% OG Author:   Mark Colton
% Edits: Spencer Jensen
% Date:     2/25/21

function ydot = mid_1_function(t,y)

% Access the global variables (model parameters) defined in the main
% function
global mu g m k r0;

% Extract the current states from the current y vector
y1 = y(1);
y2 = y(2);
theta_dot = 1.2*6.5*cos(6.5*t);
theta_ddot = -1.2*6.5*6.5*sin(6.5*t);

% Find the state derivatives
y1dot = y2;
y2dot = theta_dot^2*y1-k/m*(y1-r0)- ...
    mu*sqrt((theta_ddot*y1+2*theta_dot*y2)^2+g^2)*sign(y2);

% Reassemble the state derivatives into a single vector, ydot, to pass
% back to 'ode45'
ydot = [y1dot; y2dot];
