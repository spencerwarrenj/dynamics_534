%HW_4
%% Test Case 1:C.1.1 A quadratic Function

clear; clc; close;

x0 = [1;2]; %initial guess
tau = 1e-6; %optimality tolerance
[xstar,fstar] = minime(@objfun,x0,tau) %run with QN optimization

%run with fminunc
options = optimoptions(@fminunc,'display','iter-detailed','MaxIterations',10000,'MaxFunEvals',100000);
[xstar,fstar] = fminunc(@objfun,x0,options);


function [h,dh] = objfun(x) %objective function returns the value and derivative
    dh = [0;0]; %initialize gradient
    beta = 3/2; %set beta value (skewness)
    h = x(1)^2+x(2)^2+-beta*x(1)*x(2); %function
    dh(1) = 2*x(1)-beta*x(2); %partial with respect to x1
    dh(2) = 2*x(2)-beta*x(1); %partial with respect to x2
end