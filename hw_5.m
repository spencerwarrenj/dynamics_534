%hw_5 %newtonian dynamics

%% Problem 5.1
% syms m
m = 1;
a = 7;
b = 4;
c = 5;

d = [-a/2;-b/2;-c/2];
I_com = [m/12*(b^2+c^2),0,0;
     0,m/12*(a^2+c^2),0;
     0,0,m/12*(a^2+b^2)];
I_o_xyz = I_com+m*(d'*d*eye(3)-d*d');

theta_1 = atan2(c,a);
theta_2 = atan2(b,sqrt(a^2+c^2));
R = mrz(theta_2)*mry(-theta_1);

I_o_xyz_ddp = R*I_o_xyz*R'

%% Problem 5.2
clearvars -except I_o_xyz;clc;close;
m=.5;
L=.72;
I_1_com = [0,0,0;
           0,m*L^2/12,0;
           0,0,m*L^2/12];
I_2_com = [m*L^2/12,0,0;
           0,m*L^2/12,0;
           0,0,0];
I_3_com = I_1_com;
I_4_com = [m*L^2/12,0,0;
           0,0,0;
           0,0,m*L^2/12];
I_5_com = I_2_com;
I_6_com = I_4_com;
d1 = [-L/2;0;-L];
d2 = [-L;0;-L/2];
d3 = [-L/2;0;0];
d4 = [0;-L/2;0];
d5 = [0;-L;L/2];
d6 = [0;-L/2;L];
I_1_o = I_1_com+(d1'*d1*eye(3)-d1*d1')*m;
I_2_o = I_2_com+(d2'*d2*eye(3)-d2*d2')*m;
I_3_o = I_3_com+(d3'*d3*eye(3)-d3*d3')*m;
I_4_o = I_4_com+(d4'*d4*eye(3)-d4*d4')*m;
I_5_o = I_5_com+(d5'*d5*eye(3)-d5*d5')*m;
I_6_o = I_6_com+(d6'*d6*eye(3)-d6*d6')*m;

I_o = I_1_o+I_2_o+I_3_o+I_4_o+I_5_o+I_6_o

%% Problem 5.4
clearvars -except I_o_xyz;clc;close;
[Rt,I] = eig(I_o_xyz);
R = Rt';
R(1,:) = -1*R(1,:);
I
R
det(R)

%% Problem 5.5

I = [11.00, -2.50, -1.00;
	-2.50, 4.17, -3.33;
	-1.00, -3.33, 9.83];

[Rt,I_p] = eig(I);
R = Rt'
det(R)
R(1,:) = -1*R(1,:);
I_p
R
det(R)

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