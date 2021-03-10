%uncon
function [xstar,fstar] = minime(obj,x0,con)
%     quadratic penalty method
    
    function f_new = new_objfun(x)
        f = obj(x);
        C = con(x);
        g = max(0,C);
        f_new = f+mu_i/2*sum(g.^2);
    end

    tau = 10e-10;
    tau_feas = 10e-8;
    mu_i = 1;
    rho = 1.2;
    xstar = x0;
    x = x0;
    xstar_prev= x0-x0;
    while max(abs(xstar-xstar_prev))>tau
        
        xstar_prev = xstar;
        %minimize F
        options = optimoptions(@fminunc,'display','iter-detailed','MaxIterations',10000,'MaxFunEvals',100000);
        [xstar,fstar] = fminunc(@new_objfun,xstar,options);
        mu_i = rho*mu_i;
    end
end