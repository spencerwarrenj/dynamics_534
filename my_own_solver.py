#my_own_solver
#Author: Spencer Jensen

from grub import grub_objfun
import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint

def minime(obj,x0,con):
    #quadratic penalty method
    
    def new_objfun(x):
        nonlocal mu_e, mu_i
        f = obj(x)
        C = con(x)
        h = C[0]
        g1 = max(0,C[1])
        g2 = max(0,C[2])
        f_new = f+mu_e/2*h**2+mu_i/2*(g1**2+g2**2)
        return f_new

    tau = 10e-6
    mu_e = 1
    mu_i = 1
    rho = 1.5
    k = 0
    xstar = x0
    x = x0
    xstar_prev= x0-x0
    # options  = {'disp': True,'maxiter': 10000, 'verbose': 2}
    options  = {'disp': True}
    while np.max(np.abs(xstar-xstar_prev))>tau:
        #minimize F
        xstar_prev = xstar
        
        # lb = np.array([1e-6,0.,0.])
        # ub = np.array([10.,1.,1.])
        # bounds=Bounds(lb,ub,keep_feasible=True)
        res = minimize(new_objfun, xstar, jac='3-point' ,options=options, method='TNC') #SLSQP
        print("x = ",res.x)
        print(res.success)
        xstar = res.x
        mu_e = rho*mu_e
        mu_i = rho*mu_i
        # from IPython import embed; embed()


    fstar = obj(xstar)


    return xstar, fstar


def runoptimization(params):

    objhist = []

    joint = grub_objfun(10)

    def objcon(x):
        nonlocal objhist
        l = x[0]
        d_seg = np.array([x[1],x[2]])
        # d_seg = np.array([1])
        # print("l: ",l)
        # print("d_seg: ",d_seg)

        f = joint.get_f(d_seg,l)

        # print("f: ",f)
        g = np.zeros(3)
        g[0] = x[1]+x[2]-1
        g[1] = x[0]-1
        g[2] = -x[0]
        objhist.append(f)
        return f,g

    #------------- don't change -----------

    xlast = []
    flast = []
    glast = []

    def obj(x):
        nonlocal xlast, flast, glast
        if not np.array_equal(x,xlast):
            flast, glast = objcon(x)
            xlast = x
        return flast

    def con(x):
        nonlocal xlast, flast, glast
        if not np.array_equal(x,xlast):
            flast, glast = objcon(x)
            xlast = x
        return glast

    # ---------------------------------------

    x0 = np.array([.9,.5,.5])
    lb = np.array([1e-6,0.,0.])
    ub = np.array([10.,1.,1.])
    lg = np.array([0.,-np.inf,-np.inf])
    ug = np.array([0.,0.,0.])

    # constraints = {'type': 'ineq', 'fun': con}
    # constraints = NonlinearConstraint(con,lg,ug,jac='3-point')
    # options  = {'disp': True}
    # options  = {'disp': True,'maxiter': 10000, 'verbose': 2}
    # bounds=Bounds(lb,ub,keep_feasible=True)

#############################################
    xstar,fstar = minime(obj, x0, con) #SLSQP
#############################################

    return xstar, fstar, objhist

if __name__ == "__main__":
    params = []
    xstar,fstar, objhist = runoptimization(params)

    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogy(objhist)
    plt.show()
    plt.xlabel('Function_calls')
    plt.ylabel('Objective Function')
    
    # from IPython import embed; embed()
