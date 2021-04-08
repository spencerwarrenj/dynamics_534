%projectC
clear
%get values
get_general

%% TaskC.a


%% Task C.b
clear;
clc;

% Time (s)
dt = 0.01;      % Step size
tf = 30;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 0;
wy0 = 0;
wz0 = 0;
psi0 = 0;
theta0 = 0;
phi0 = 0;

% Torques (N-m). You can change these torque vectors to validate your
% model and try out the different cases required in the project. When
% you send me your function, I will try out some torques to see if
% your model accurately predicts the response.
Mx = 176*cos(.2*t); 
% Mx(1:50) = -176;
My = 54*ones(size(t));  %zeros(size(t)); %
Mz = 98*sin(0.3*t); %zeros(size(t)); %

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=jensen(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

max(wx)
max(wy)
max(phi)

% Plot results
subplot(2,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z');
subplot(2,1,2);
plot(t,psi,t,theta,t,phi);
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi');

%% Task C.c moments
%find Moment to apply for barbecue mode
I = [40823.1019533526,1537.80758482821,-3179.29746413299;1537.80758482821,90593.5331952687,128.577121406529;-3179.29746413299,128.577121406529,98742.9021241113];
omega = [1;0;0]*pi/180;
tau = cross(omega,I*omega)

%% Task C.c case1
%simulate barbecue mode
clearvars -except tau;
clc;

% Time (s)
dt = 0.01;      % Step size
tf = 1000;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 1;
wy0 = 0;
wz0 = 0;
psi0 = 0;
theta0 = 0;
phi0 = 0;

% Torques (N-m). You can change these torque vectors to validate your
% model and try out the different cases required in the project. When
% you send me your function, I will try out some torques to see if
% your model accurately predicts the response.
Mx = tau(1)*ones(size(t));
My = tau(2)*ones(size(t));
Mz = tau(3)*ones(size(t));

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=jensen(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

% Plot results
subplot(2,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z');
subplot(2,1,2);
plot(t,psi,t,theta,t,phi);
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi');

%% Task C.c case2
%simulate no torques
clearvars -except tau;
clc;

% Time (s)
dt = 0.01;      % Step size
tf = 1000;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 1;
wy0 = 0;
wz0 = 0;
psi0 = 0;
theta0 = 0;
phi0 = 0;

% Torques (N-m). You can change these torque vectors to validate your
% model and try out the different cases required in the project. When
% you send me your function, I will try out some torques to see if
% your model accurately predicts the response.
Mx = zeros(size(t));
My = zeros(size(t));
Mz = zeros(size(t));

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi]=jensen(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

% Plot results
subplot(2,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z');
subplot(2,1,2);
plot(t,psi,t,theta,t,phi);
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi');

