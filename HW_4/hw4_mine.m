function [xopt, fopt, iter, opt] = hw_3_comp_step()
clear; clc; close;
x0 = ones(1,10); %initial guess
lb = ones(1,10).*.1; %lower bound
ub = []; %not used
y_stress = ones(10,1).*25*10^3; %array of yield stress
y_stress(9) = 75*10^3; %change the 9th yield stress
f_count = 0;
function [f, g, h, df, dg, dh] = objcon(x)
    f_count = f_count+1;
    [mass,stress] = truss(x); %evaluate function
    f = mass; %objective function
    f = f./100;
    g = zeros(20,1);
    g(1:10,1) = abs(stress)-y_stress; %calculate the inequality constraint
    g(1:10,1) = g(1:10,1)./y_stress;
    g(11:20,1) = -x+.1;
    h = []; %unused
    
    step = 10e-30; %step length for finite diff
    %get df
    df = zeros(length(x),1);
    for i = 1 : length(x)
        newx = x; %initialize new x
        newx(i) = x(i)+1i*step; %take step in ith direction
        [dmass,~] = truss(newx); %evaluate new x
        df(i) = imag(dmass)/step; %calculate slope
    end
    
    %get dg
    nx = [10]; %number of design variables
    ng = [10]; %number of contraints
    dg = zeros(ng,nx); %initialize dg
    for i = 1 : length(x)
        newx = x; %set x
        newx(i) = x(i)+1i*step; %take a step in the ith design variable
        [~,stress] = truss(newx); %evaluate stress
        dg(:,i) = imag(stress)/step.*(sign(real(stress))); %difference
    end
    dg = dg'; %inequality constraints
    
    dh = []; %equality constraints
end

% -----don't change -------------------------------------
                       
% ------- shared variables -----------
xlast = [];  % last design variables
flast = [];  % last objective
glast = [];  % last nonlinear inequality constraint
hlast = [];  % last nonlinear equality constraint
dflast = [];  % last derivatives
dglast = [];
dhlast = [];
% --------------------------------------

% ------ separate obj/con  --------
function [f] = obj(x)

    % check if computation is necessary
    if ~isequal(x, xlast)
        [flast, glast, hlast, dflast, dglast, dhlast] = objcon(x);
        xlast = x;
    end

    f = flast;
    df = dflast;
end

function [g] = con(x)

    % check if computation is necessary
    if ~isequal(x, xlast)
        [flast, glast, hlast, dflast, dglast, dhlast] = objcon(x);
        xlast = x;
    end

    % set constraints
    g = glast;
    h = hlast;
    dg = dglast;
    dh = dhlast;
end
% ------------------------------------------------

% call fmincon
[xopt, fopt] = minime(@obj, x0, @con);

f_count
end

