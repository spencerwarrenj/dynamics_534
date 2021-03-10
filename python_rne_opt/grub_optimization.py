from grub import grub_objfun
import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint


def runoptimization(params):

    objhist = []

    joint = grub_objfun(100)

    def objcon(x):
        nonlocal objhist
        l = x[0]
        d_seg = np.array([x[1],x[2]])
        # d_seg = np.array([1])
        # print("l: ",l)
        # print("d_seg: ",d_seg)

        f = joint.get_f(d_seg,l)

        # print("f: ",f)
        g = np.zeros(2)
        g[0] = x[1]+x[2]-1
        g[1] = x[0]-1
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
    # x0 = np.array([.1,.4,.6])
    lb = np.array([1e-6,0.,0.])
    ub = np.array([10.,1.,1.])
    lg = np.array([0.,-np.inf])
    ug = np.array([0.,0.])

    # constraints = {'type': 'ineq', 'fun': con}
    constraints = NonlinearConstraint(con,lg,ug,jac='3-point')
    # options  = {'disp': True}
    options  = {'disp': True,'maxiter': 10000, 'verbose': 2}
    bounds=Bounds(lb,ub,keep_feasible=True)

    res = minimize(obj, x0, jac='3-point', constraints=constraints, options=options, bounds=bounds, method='trust-constr') #SLSQP
    print("x = ",res.x)
    print("f = ",res)
    print(res.success)

    return res.x, res.fun, objhist

if __name__ == "__main__":
    params = []
    xstar,fstar, objhist = runoptimization(params)

    print("xstar: ",xstar)
    print("fstar: ",fstar)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogy(objhist)
    
    plt.xlabel('Function_calls')
    plt.ylabel('Objective Function')
    plt.show()
    
    # from IPython import embed; embed()