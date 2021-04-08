from python_rne_opt import grub_rne as dyn
import numpy as np
from skopt.space import Space
from skopt.sampler import Lhs

class grub_objfun:
    def __init__(self,positions):
        self.positions = positions #number of random configurations to try

        np.random.seed(3)

        n_samples = self.positions
        lhs = Lhs(criterion="maximin", iterations=10000)
        space = Space([(-1.57, 1.57),(-1.57,1.57),(-2., 2.),(-2.,2.),(-2., 2.),(-2.,2.)])
        x = lhs.generate(space.dimensions, n_samples)
        x = np.array(x)
        self.uv_mat = x[:,:2]
        self.uvdot_mat = x[:,2:4]
        self.tau_mat = x[:,4:]
        print("uvs: ",self.uv_mat)

        self.true_xDot = np.zeros((self.positions,2))
        for i in range(self.positions):
            #set uv, uvdot, tau        
            uv = self.uv_mat[i,0:].flatten()
            uvdot = self.uvdot_mat[i,0:].flatten()
            tau = self.tau_mat[i,0:].flatten()

            #get "true answer"
            d_seg_act = np.ones(1000)
            d_seg_act = d_seg_act*1/1000
            l_act = 1
            grub = dyn.rne_model("joint4",d_seg_act,l_act)
            self.true_xDot[i,:] = grub.get_xdot(uv,uvdot,tau)
            # print("True qdd: ",self.true_xDot[i,:])
        
    def get_f(self,d_seg,l):
        f = np.zeros(self.positions)
        # print(f)
        for i in range(self.positions):
            #set uv, uvdot, tau
            uv = self.uv_mat[i,0:]
            uvdot = self.uvdot_mat[i,0:]
            tau = self.tau_mat[i,0:]

            #find difference between true and approx.
            grub = dyn.rne_model("joint4",d_seg,l)
            xDot = grub.get_xdot(uv,uvdot,tau)
            # print("xDot: ",xDot)
            # print("f: ",np.linalg.norm(xDot-true_xDot))
            # from IPython import embed; embed()
            f[i] = np.linalg.norm((xDot-self.true_xDot[i,:]))
            # f[i] = np.max(np.abs(xDot-self.true_xDot[i,:])/self.true_xDot[i,:])
            # print("f: ",f[i])
        mean_f = np.max(f)
        return mean_f

##Testing grub function
if __name__ == "__main__":
    d_seg_act = np.ones(1) 
    d_seg_act = d_seg_act*1/1
    l_act = 1
    grub = dyn.rne_model("joint4",d_seg_act,l_act)
    uv = np.array([1,1])
    uvdot = np.array([1,1])
    tau = np.array([0,0]) 
    xdot = grub.get_xdot(uv,uvdot,tau)
    print("xdot: ",xdot)
