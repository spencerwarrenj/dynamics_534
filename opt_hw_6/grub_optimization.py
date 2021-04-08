from grub import grub_objfun
import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint
from scipy.optimize import differential_evolution,dual_annealing
import matplotlib.pyplot as plt
import time

def runoptimization(num_seg,method):

    objhist = []
    # plt.ion()
    # plt.show()

    joint = grub_objfun(100)

    def objcon(x):
        nonlocal objhist,method
        l = x[0]
        d_seg = np.array([x[1:]])
        # d_seg = np.array([1])
        # print("l: ",l)
        # print("d_seg: ",d_seg)

        f = joint.get_f(d_seg.flatten(),l)

        #for nelder-mead
        if method == "nm":
            if x[0] < 0:
                f = f+10000*x[0]**2
            if x[0] > 1:
                f = f+10000*x[0]**2
            mu = 10000000
            g_pm = np.sum(x[1:])-1
            f = f+mu/2*g_pm**2


        # print("f: ",f)
        g = np.zeros(1)
        g[0] = np.sum(x[1:])-1
        objhist.append(f)

        # plt.semilogy(objhist)
        # plt.draw()
        # plt.pause(0.001)
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
    def seeing(xk,convergence=1):
        plt.semilogy(objhist,'b')
        plt.draw()
        plt.xlabel('Function_calls')
        plt.ylabel('Objective Function')
        plt.pause(0.001)


    # ---------------------------------------

    x0 = np.array([.9])
    x0_segs = np.ones(num_seg)
    x0 = np.append(x0,x0_segs*1/num_seg)
    # x0 = np.array([.9,.5,.5])
    # x0 = np.array([.1,.4,.6])
    lb = np.array([.01])
    ub = np.array([1.])
    lb = np.append(lb,np.zeros(num_seg))
    ub = np.append(ub,np.ones(num_seg))
    # from IPython import embed; embed()
    # lb = np.array([.01,0.,0.])
    # ub = np.array([1.,1.,1.])
    lg = np.array([0.])
    ug = np.array([0.])

    # constraints = {'type': 'ineq', 'fun': con}
    constraints = NonlinearConstraint(con,lg,ug,jac='3-point')
    # options  = {'disp': True}
    options  = {'disp': True,'maxiter': 10000, 'verbose': 2}
    bounds=Bounds(lb,ub,keep_feasible=True)

    nlc = NonlinearConstraint(con,lg,ug)

    if method == "tr":
        res1 = minimize(obj, x0, jac='3-point', constraints=constraints, options=options, bounds=bounds, method='trust-constr') #SLSQP
        x0 = np.array([.9])
        left = 1
        num = 1
        for i in range(num_seg):
            if i == num_seg-1:
                x0 = np.append(x0,left)
            else:
                new_num = num*.5
                x0 = np.append(x0,new_num)
                left = 1-sum(x0[1:])
                num = new_num
        res2 = minimize(obj, x0, jac='3-point', constraints=constraints, options=options, bounds=bounds, method='trust-constr') #SLSQP
        if res1.fun<res2.fun:
            res = res1
        else:
            res = res2
    if method == "nm":
        res1 = minimize(obj, x0, jac='3-point', constraints=constraints, options=options, bounds=bounds, method='Nelder-Mead') #SLSQP
        x0 = np.array([.9])
        left = 1
        num = 1
        for i in range(num_seg):
            if i == num_seg-1:
                x0 = np.append(x0,left)
            else:
                new_num = num*.5
                x0 = np.append(x0,new_num)
                left = 1-sum(x0[1:])
                num = new_num
        # from IPython import embed; embed()
        res2 = minimize(obj, x0, jac='3-point', constraints=constraints, options=options, bounds=bounds, method='Nelder-Mead') #SLSQP
        if res1.fun<res2.fun:
            res = res1
        else:
            res = res2
    if method == "ga":
        res = differential_evolution(obj, bounds, callback=seeing, constraints=nlc,seed=1,polish=False)
    # if method == "da":
    #     res = dual_annealing(obj,bounds=bounds)
    print("x = ",res.x)
    print("f = ",res)
    print(res.success)

    # plt.savefig("evolutionary.png") #, dpi=None, facecolor='w', edgecolor='w',
    #     orientation='portrait', papertype=None, format=None,
    #     transparent=False, bbox_inches=None, pad_inches=0.1,
    #     frameon=None, metadata=None)

    return res.x, res.fun, objhist

if __name__ == "__main__":
    xstar_ga = [None] * 10
    fstar_ga = np.zeros(10)
    for i in range(10):
        start = time.perf_counter()
        num_seg = i+1
        method = 'ga'
        xstar_ga[i],fstar_ga[i], objhist = runoptimization(num_seg,method)
        print("xstar: ",xstar_ga)
        print("fstar: ",fstar_ga)
        end = time.perf_counter()  
        time_ga = end-start
        print("took " + str(end-start) + " Seconds")
    
    np.save('optimal_ga.npy', (fstar_ga,xstar_ga), allow_pickle=True)
    # np.load('optimal_ga.npy', allow_pickle=True)

    xstar_nm = [None] * 10
    fstar_nm = np.zeros(10)
    for i in range(10):
        start = time.perf_counter()
        num_seg = i+1
        method = 'nm'
        xstar_nm[i],fstar_nm[i], objhist = runoptimization(num_seg,method)
        print("xstar: ",xstar_nm)
        print("fstar: ",fstar_nm)
        end = time.perf_counter()  
        time_nm = end-start
        print("took " + str(end-start) + " Seconds")
    
    np.save('optimal_nm.npy', (fstar_nm,xstar_nm), allow_pickle=True)
    # np.load('optimal_nm.npy', allow_pickle=True)

    xstar_tr = [None] * 10
    fstar_tr = np.zeros(10)
    for i in range(10):
        start = time.perf_counter()
        num_seg = i+1
        method = 'tr'
        xstar_tr[i],fstar_tr[i], objhist = runoptimization(num_seg,method)
        print("xstar: ",xstar_tr)
        print("fstar: ",fstar_tr)
        end = time.perf_counter()  
        time_tr = end-start
        print("took " + str(end-start) + " Seconds")
    np.save('optimal_tr.npy', (fstar_tr,xstar_tr), allow_pickle=True)
    # np.load('optimal_tr.npy', allow_pickle=True)

    from IPython import embed; embed()
    
    # plt.figure()
    # plt.semilogy(objhist)
    
    # plt.xlabel('Function_calls')
    # plt.ylabel('Objective Function')
    # # plt.savefig("nelder_mead.png")
    # plt.show()
    
    
    # from IPython import embed; embed()