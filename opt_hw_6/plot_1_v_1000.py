#plotting RNE against actual data

from python_rne_opt import grub_rne as dyn
import numpy as np
import scipy.io
import scipy as sp
import scipy.integrate as spint
import matplotlib.pyplot as pplot
import matplotlib

def plot_stuff(num_seg):

    tau_hist = np.array([])
    uvdot_hist = np.array([])

    def closest(lst, K):
        lst = np.asarray(lst)
        idx = (np.abs(lst - K)).argmin()
        return idx

    def fun_sim(t,y):
        nonlocal t_vec,p0,p1,p2,p3,tau_hist,uvdot_hist,num_seg
        d = .15
        r = .16
        A = np.pi*r**2
        ind = closest(t_vec,t)
        # tau = np.array([d*A*(p0[ind]-p1[ind]),d*A*(p2[ind]-p3[ind])])
        tau = np.array([0.,0.])
        d_seg_act = np.ones(num_seg) 
        d_seg_act = d_seg_act*1/num_seg
        l_act = 1
        grub = dyn.rne_model("joint4",d_seg_act,l_act)
        uv = y[:2]
        uvdot = y[2:]
        uvddot = grub.get_xdot(uv,uvdot,tau)
        # from IPython import embed; embed()
        tau_hist = np.append(tau_hist,tau)
        uvdot_hist = np.append(uvdot_hist,uvdot)
        print("t: ",t)
        # from IPython import embed; embed()
        return np.concatenate((uvdot,uvddot))

    mat = scipy.io.loadmat('FULL_SIM+ERR_hw_data_Thu_Mar_11_14-41-02_2021.mat')
    start_ind = 70
    end_ind = 130
    u = mat['x_hist'][6,start_ind:end_ind].flatten()
    v = mat['x_hist'][7,start_ind:end_ind].flatten()
    t = mat['timestamp'][:,start_ind:end_ind].flatten()
    t_vec = t
    # pplot.plot(t,u,'b-',t,v,'r-')
    

    udot = mat['x_hist'][4,start_ind:end_ind].flatten()
    vdot = mat['x_hist'][5,start_ind:end_ind].flatten()
    p0 = mat['x_hist'][0,start_ind:end_ind].flatten()
    p1 = mat['x_hist'][1,start_ind:end_ind].flatten()
    p2 = mat['x_hist'][2,start_ind:end_ind].flatten()
    p3 = mat['x_hist'][3,start_ind:end_ind].flatten()
    # from IPython import embed; embed()
    print(np.size(u))
    # uv0 = np.array([[u[0]],[v[0]]])
    # uvdot0 = np.array([[udot[0]],[vdot[0]]])
    uv0 = np.array([[.01],[-.01]])
    uvdot0 = np.array([[0.],[0.]])
    y0 = np.vstack((uv0,uvdot0))
    print(y0)
    t_span = (0,1)
    res = spint.solve_ivp(fun_sim, t_span, y0.flatten(), method='RK23')
    # from IPython import embed; embed()
    return res


if __name__ == "__main__":
    pplot.rcParams.update({'font.size': 14})

    num_seg = 1
    res = plot_stuff(num_seg)
    u = res.y[0,:]
    v = res.y[1,:]
    pplot.plot(res.t,u,'b-',label='1-segment')
    pplot.plot(res.t,v,'r-')
    # pplot.legend(['u 1-segment','v 1-segment'])

    # pplot.show()

    num_seg = 10
    res = plot_stuff(num_seg)
    u = res.y[0,:]
    v = res.y[1,:]
    pplot.plot(res.t,u,'b-.',label='10-segment')
    pplot.plot(res.t,v,'r-.')

    num_seg = 100
    res = plot_stuff(num_seg)
    u = res.y[0,:]
    v = res.y[1,:]
    pplot.plot(res.t,u,'b--',label='100-segment')
    pplot.plot(res.t,v,'r--')

    num_seg = 1000
    res = plot_stuff(num_seg)
    u = res.y[0,:]
    v = res.y[1,:]
    pplot.plot(res.t,u,'b:',label='1000-segment')
    pplot.plot(res.t,v,'r:')
    # pplot.legend(['u 1-segment','v 1-segment','u 10-segment','v 10-segment','u 100-segment','v 100-segment','u 1000-segment','v 1000-segment'])
    pplot.legend()
    pplot.xlabel('Time (s)')

    pplot.ylabel('Rotation (rad)')
    


    pplot.show()
