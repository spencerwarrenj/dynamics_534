function [xopt, fopt, iter, opt] = hw_3_AD()
clear; clc; close;
x0 = ones(1,10); %initial guess
lb = ones(1,10).*.1; %lower bound
ub = []; %not used
y_stress = ones(10,1).*25*10^3; %array of yield stress
y_stress(9) = 75*10^3; %change the 9th yield stress

function [f, g, h, df, dg, dh] = objcon(x)
    
    [mass,stress] = truss(x); %evaluate function
    f = mass; %objective function
    g = abs(stress)-y_stress; %calculate the inequality constraint
%     xnew = amatinit(x);
%     [mass,stress] = truss(xnew); %evaluate function
%     dstress = abs(stress);
    
    
    
    h = []; %unused
    
%     df = ajac(mass,0);
    df = [];
%     dg = ajac(dstress,0);
%     dg = dg';
    dg = [];
    dh = []; %equality constraints
end

% -----don't change -------------------------------------

options = optimoptions(@fmincon, ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', 1e5, ...
        'SpecifyObjectiveGradient', false, ...
        'SpecifyConstraintGradient', false, ...
        'OutputFcn', @outfun); %, ...
%         'CheckGradients', true);

iter = [];
opt = [];
func_calls = [];

function stop = outfun(x,ov,state)
    iter = [iter; ov.iteration];
    opt = [opt; ov.firstorderopt];
    func_calls = [func_calls; ov.funccount];
    stop = false;
end

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
function [f, df] = obj(x)

    % check if computation is necessary
    if ~isequal(x, xlast)
        [flast, glast, hlast, dflast, dglast, dhlast] = objcon(x);
        xlast = x;
    end

    f = flast;
    df = dflast;
end

function [g, h, dg, dh] = con(x)

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

% linear constraints not used
A = [];
b = [];
Aeq = [];
beq = [];
% call fmincon
[xopt, fopt] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);

figure(1)
semilogy(func_calls,opt)
xlabel('Number of Function Calls')
ylabel('First-order Optimality')
set(gca,'Fontsize',12);
saveas(figure(1),'convergence_plot.jpg')

end