% This is the function that 'ode45' calls to get the function derivatives
% at each time step.  'ode45' passes the current time and states (x and v,
% which are contained in y), and the function returns the current state
% derivatives (xdot and vdot, which are contained in ydot).
%
% OG Author:   Mark Colton
% Edits: Spencer Jensen
% Date:     2/18/21

function ydot = hw3_5_function(t,y)

% Access the global variables (model parameters) defined in the main
% function
global m k R;

% Extract the current states from the current y vector
y1 = y(1);
y2 = y(2);
w = 5*sin(t);
w_dot = 5*cos(t);

% Find the state derivatives
y1dot = y2;
y2dot = -k*y1/m+R*w_dot+y1*w^2;

% Reassemble the state derivatives into a single vector, ydot, to pass
% back to 'ode45'
ydot = [y1dot; y2dot];


% ALTERNATIVE, ALL-IN-ONE APPROACH IN WHICH I DON'T EXPLICITLY SEPARATE
% OUT THE STATES.  IT IS IMPLICIT THAT y(1) = x AND y(2) = v.
% ydot = [y(2); 1/m/(y(1)^2+1)*(m*(omega^2-g)*y(1)-m*y(1)*y(2)^2-k/2*y(1)^3)];
