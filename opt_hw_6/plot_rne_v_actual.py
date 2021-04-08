#plotting RNE against actual data

from python_rne_opt import grub_rne as dyn
import numpy as np
import scipy.io
import scipy as sp
import scipy.integrate as spint
import matplotlib.pyplot as pplot

def plot_stuff():

    tau_hist = np.array([])
    uvdot_hist = np.array([])

    def closest(lst, K):
        lst = np.asarray(lst)
        idx = (np.abs(lst - K)).argmin()
        return idx

    def fun_sim(t,y):
        nonlocal t_vec,p0,p1,p2,p3,tau_hist,uvdot_hist
        d = .15
        r = .16
        A = np.pi*r**2
        ind = closest(t_vec,t)
        tau = np.array([d*A*(p0[ind]-p1[ind]),d*A*(p2[ind]-p3[ind])])
        d_seg_act = np.ones(1) 
        d_seg_act = d_seg_act*1/1
        l_act = 1
        grub = dyn.rne_model("joint4",d_seg_act,l_act)
        uv = y[:2]
        uvdot = y[2:]
        uvddot = grub.get_xdot(uv,uvdot,tau)
        # from IPython import embed; embed()
        tau_hist = np.append(tau_hist,tau)
        uvdot_hist = np.append(uvdot_hist,uvdot)
        # print("tau: ",tau)
        # from IPython import embed; embed()
        return np.concatenate((uvdot,uvddot))

    mat = scipy.io.loadmat('FULL_SIM+ERR_hw_data_Thu_Mar_11_14-41-02_2021.mat')
    start_ind = 10
    end_ind = 8000
    u = mat['x_hist'][6,start_ind:end_ind].flatten()
    v = mat['x_hist'][7,start_ind:end_ind].flatten()
    t = mat['timestamp'][:,start_ind:end_ind].flatten()
    t_vec = t
    pplot.plot(t,u,'b-',t,v,'r-')
    pplot.legend(['u','v'])

    udot = mat['x_hist'][4,start_ind:end_ind].flatten()
    vdot = mat['x_hist'][5,start_ind:end_ind].flatten()
    p0 = mat['x_hist'][0,start_ind:end_ind].flatten()
    p1 = mat['x_hist'][1,start_ind:end_ind].flatten()
    p2 = mat['x_hist'][2,start_ind:end_ind].flatten()
    p3 = mat['x_hist'][3,start_ind:end_ind].flatten()
    # from IPython import embed; embed()
    print(np.size(u))
    uv0 = np.array([[u[0]],[v[0]]])
    uvdot0 = np.array([[udot[0]],[vdot[0]]])
    y0 = np.vstack((uv0,uvdot0))
    print(y0)
    t_span = (t[0],t[-1])
    res = spint.solve_ivp(fun_sim, t_span, y0.flatten(), method='RK45')
    # from IPython import embed; embed()
    u = res.y[0,:]
    v = res.y[1,:]
    pplot.plot(res.t,u,'b:',res.t,v,'r:')
    
    pplot.show()

    print("Max tau: ",max(tau_hist))
    print("Min tau: ",min(tau_hist))
    pplot.figure()
    pplot.hist(uvdot_hist)
    pplot.title('uvdot')
    pplot.figure()
    pplot.hist(tau_hist)
    pplot.title('tau')

    pplot.show()


if __name__ == "__main__":
    plot_stuff()
# scipy.integrate.odeint(func, y0, t, args=()
# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', args=None)