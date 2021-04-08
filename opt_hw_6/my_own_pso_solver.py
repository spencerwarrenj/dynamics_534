from grub import grub_objfun
import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt




def runoptimization(params):

    objhist = []
    # plt.ion()
    # plt.show()

    joint = grub_objfun(100)

    def minime_pso(obj):
        #particle swarm optimization
        nonlocal objhist
        #initialize population
        pop_size = 100
        design_vars = 2
        x = np.random.rand(pop_size,design_vars)
        x_best_not_changed = 0
        omega = .9
        c1 = np.random.rand(pop_size,1)*2
        c2 = np.random.rand(pop_size,1)*2
        delta_x = np.random.rand(pop_size,design_vars)*2
        x_i_best = x
        
        # find best member
        f_best = 100000
        for i in range(pop_size):
            f = obj(x[i,:])
            if f < f_best:
                f_best = f
                x_best = x[i,:]

        #set convergence criteria
        delta_x_max = 1
        while x_best_not_changed < 100: 
            # from IPython import embed; embed()
            delta_x = omega*delta_x+c1*(x_i_best-x)+c2*(x_best-x)
            # from IPython import embed; embed()
            delta_x_max = .9*delta_x_max
            delta_x[delta_x>delta_x_max] = delta_x_max
            x = x + delta_x
            x[x<0]=0
            #find population best
            f_new_best = 1000000
            for i in range(pop_size):
                f = obj(x[i,:])
                if f < f_new_best:
                    f_new_best = f
                    x_new_best = x[i,:]
            #store population best if better than best
            if f_new_best < f_best:
                f_best = f_new_best
                x_best = x_new_best
                print("delta_x_max: ",delta_x_max)
                # from IPython import embed; embed()
                plt.semilogy(objhist,'b')
                plt.draw()
                plt.xlabel('Function_calls')
                plt.ylabel('Objective Function')
                plt.pause(0.001)
            else:
                x_best_not_changed = x_best_not_changed+1

            
        xstar = x_best
        #return best objective function
        fstar = obj(xstar)


        return xstar, fstar


    def objcon(x):
        nonlocal objhist
        l = x[0]
        d_seg = np.array([x[1],1-x[1]])

        f = joint.get_f(d_seg,l)
        objhist.append(f)
        return f

    #------------- don't change -----------

    xlast = []
    flast = []
    glast = []

    def obj(x):
        nonlocal xlast, flast, glast
        if not np.array_equal(x,xlast):
            flast = objcon(x)
            xlast = x
        return flast

    xstar,fstar = minime_pso(obj)

    return xstar, fstar, objhist

if __name__ == "__main__":
    params = []
    xstar,fstar, objhist = runoptimization(params)

    print("xstar: ",xstar)
    print("fstar: ",fstar)
    
    plt.figure()
    plt.semilogy(objhist)
    
    plt.xlabel('Function_calls')
    plt.ylabel('Objective Function')
    # plt.savefig("nelder_mead.png")
    plt.show()
    
    
    # from IPython import embed; embed()