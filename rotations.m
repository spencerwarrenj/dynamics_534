%HW 1
%% problem 1
clear
%part a
syms phi theta psi
R = mrz(psi)*mrx(theta)*mrz(phi)
%part b
phi = 30*pi/180;
theta = 45*pi/180;
psi = -60*pi/180;
R = double(simplify(subs(R)))
%part c
r_xyz = [-5;3;0];
r_XYZ = R'*r_xyz

%% problem 2 (3.3 from book)
clear
r_airplane = [0;.5;-2];
R = mry(10*pi/180)*mrx(20*pi/180)*mrz(40*pi/180);
r_earth = R*r_airplane

%% problem 3 (3.4 from book)
clear;
e_ba = [-50;20;0];
e_ca = [-50;0;40];
e_bah = e_ba/norm(e_ba);
e_cah = e_ca/norm(e_ca);
k_hat = cross(e_bah,e_cah);
i_hat = e_bah;
j_hat = cross(k_hat,i_hat);
R = [i_hat';j_hat';k_hat']
r_XYZ = [-50;0;0];
r_xyz = R*r_XYZ

%% problem 4
clear;
syms theta
R_rot = mry(theta)*mrz(40*pi/180)
theta_r = double(solve(R_rot(1,1)==cos(49.36*pi/180),theta))
theta_r*180/pi

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