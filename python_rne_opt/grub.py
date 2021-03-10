from python_rne_opt import grub_rne as dyn
import numpy as np

class grub_objfun:
    def __init__(self,positions):
        self.positions = positions #number of random configurations to try

        np.random.seed(3)
        self.uv_mat = np.random.rand(self.positions,2)*2-1
        self.uvdot_mat = np.random.rand(self.positions,2)*2-1
        self.tau_mat = np.random.rand(self.positions,2)*2-1
        print("uvs: ",self.uv_mat)

        self.true_xDot = np.zeros((self.positions,2))
        for i in range(self.positions):
            #set uv, uvdot, tau        
            uv = self.uv_mat[i,0:].flatten()
            uvdot = self.uvdot_mat[i,0:].flatten()
            tau = self.tau_mat[i,0:].flatten()

            #get "true answer"
            d_seg_act = np.ones(100)
            d_seg_act = d_seg_act*1/100
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
            f[i] = np.linalg.norm(xDot-self.true_xDot[i,:])
            # print("f: ",f[i])
        mean_f = np.max(f)
        return mean_f