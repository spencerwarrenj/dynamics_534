%hw4
%% problem 2
clear; clc; close;
syms T theta Na Nb theta_dot theta_ddot L m g mu
% eq1 = Na==T*cos(theta)
eq2 = T*cos(theta)-m*g==m*(-L*theta_ddot*sin(theta)-L*theta_dot^2*cos(theta))
eq3 = T*sin(theta)-Nb*mu == m*(L*theta_ddot*cos(theta)-L*theta_dot^2*sin(theta))
eq4 = Nb-T*cos(theta)-m*g==0
Nb_act = solve(eq4,Nb)
eq5 = subs(eq3,Nb,Nb_act)
T_act = solve(eq2,T)
my_sol = simplify(subs(eq5,T,T_act))


solution = theta_ddot*(cos(theta)*(1+tan(theta)^2)-mu*sin(theta))-theta_dot^2*mu*cos(theta)-g/L*tan(theta)+g*mu*g/L

simplify(my_sol-solution)