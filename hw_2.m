%HW_2

%% HW_2.1
psi = 60*pi/180;
psi_dot = -.3;
w = 1500*2*pi/60;
omega = [psi_dot;w*sin(psi);w*cos(psi)]
alpha = [0;w*psi_dot*cos(psi);-w*psi_dot*sin(psi)]

%% HW_2.2
w_3 = 5000*2*pi/60;
w_1 = 3;
w_1_dot = -1.8;
theta = -75;
theta_ddot = -3;
omega_dot_rel = [w_1_dot;theta_ddot;0];
omega_frame = [w_1;0;0];
omega_extra = [w_3*cos(theta);0;w_3*sin(theta)];
% alpha = omega_dot_rel+[0;-w_1*w_3*sin(theta);0];
alpha = omega_dot_rel+cross(omega_frame,omega_extra)

%% HW_2.3 (prob. 3.38 from text with modifications)
syms w1 w2 beta psi L1 L2 R
w = [-w1*sin(beta);0;w1*cos(beta)];
r = [L1*cos(beta)+L2+R*cos(psi);R*sin(psi);L1*cos(beta)];
cross_term = cross(r,w) 
