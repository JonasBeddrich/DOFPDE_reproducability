import matplotlib.pyplot
import numpy as np 
from scipy.optimize import fsolve
import matplotlib.pyplot as plt 
from scipy.special import gamma
import matplotlib
from tqdm import * 

from mpmath import hyper

class FDE02RK: 
    
    # fixed dt 
    # scalar y 

    def __init__(self, weights01, exponents01, weights12, exponents12, f, t0, u0, v0, tf, h, scheme="GL1"): 

        # fractional part 
        self.m01 = len(weights01)
        self.m12 = len(weights12)

        self.weights01 = weights01
        self.exponents01 = exponents01
        self.weights12 = weights12
        self.exponents12 = exponents12

        # right hand side 
        self.f = f 

        # initial datum 
        self.t0 = t0 
        self.u0 = u0
        self.v0 = v0 

        # algorithm 
        self.tf = tf 
        self.h = h 
        # only fixed steps right now 

        # Set Butcher Table 
        self.selectRK(scheme)

    def solve(self): 
        t_values = np.arange(self.t0, self.tf+1e-12, self.h)
        # initial condition of modes is 0 
        y_values = np.zeros((len(t_values), 2 + self.m01 + self.m12)) 
        # initial condition of u and v
        y_values[0,0] = self.u0
        y_values[0,1] = self.v0 

        self.preprocess()

        for i in tqdm(range(1, len(t_values))):
            y_values[i] = self.step(t_values[i-1], y_values[i-1])

        return t_values, y_values

    def step(self, tn, yn):

        y_next = np.empty(yn.shape)

        unc = self.solveNonlinearEquation(tn, yn)
        vnc = self.computeVNC(un=yn[0], unc=unc)
        rnc = self.computeRNC(vn=yn[1], vnc=vnc)
        
        y_next[0] = self.computeSolution(un=yn[0], vnc=vnc)
        y_next[1] = self.computeDerivative(vn = yn[1], rnc=rnc)
        y_next[2:2+self.m01] = self.computeModes01(yn=yn, vnc=vnc)
        y_next[2+self.m01:] = self.computeModes12(yn=yn, rnc=rnc) 

        return y_next

    def computeVNC(self, un, unc): 
        return np.matmul(self.Ainv, unc - un * self.vec1) / self.h

    def computeRNC(self, vn, vnc): 
        return np.matmul(self.Ainv, vnc - vn * self.vec1) / self.h 

    def computeSolution(self, un, vnc): 
        return un + self.h * np.sum(vnc * self.b)
    
    def computeDerivative(self, vn, rnc): 
        return vn + self.h * np.sum(rnc * self.b)

    def computeModes01(self, yn, vnc):
        ak_steps = np.empty((self.m01, self.RKsteps))
        for k in range(self.m01): 
            ak_steps[k] = np.matmul(self.Bk[k], - self.exponents01[k] * yn[2+k] * self.vec1 + self.weights01[k] * vnc)
        return yn[2:2+self.m01] + self.h * np.matmul(ak_steps, self.b)

    def computeModes12(self, yn, rnc): 
        bj_steps = np.empty((self.m12, self.RKsteps))
        for j in range(self.m12): 
            bj_steps[j] = np.matmul(self.Cj[j], - self.exponents12[j] * yn[2+self.m01+j] * self.vec1 + self.weights12[j] * rnc)
        return yn[2+self.m01:] + self.h * np.matmul(bj_steps, self.b)

    def solveNonlinearEquation(self, tn, yn): 
        tmp = np.zeros(self.RKsteps)
        # from a_k 
        tmp += np.sum(yn[2:self.m01+2]) * self.vec1
        for k in range(self.m01): 
            tmp -= self.h * self.exponents01[k] * yn[2+k] * np.matmul(self.A, np.matmul(self.Bk[k], self.vec1))
        tmp -= yn[0] * np.matmul(self.D, self.vec1)
        
        # from b_j 
        tmp += np.sum(yn[self.m01+2:]) * self.vec1
        for j in range(self.m12): 
            tmp -= self.h * self.exponents12[j] * yn[self.m01+2+j] * np.matmul(self.A, np.matmul(self.Cj[j], self.vec1))
        tmp -= yn[0] * np.matmul(self.F, self.vec1)
        tmp -= yn[1] * np.matmul(self.E, self.vec1)     

        def nonlinearequation(unc): 
            r = tmp + 0 
            r += np.matmul(self.D + self.F, unc)
            for i in range(self.RKsteps): 
                r[i] -= self.f(tn + self.c[i] * self.h, unc[i])
            return r 
        result = fsolve(nonlinearequation, yn[0] * np.ones(self.RKsteps), xtol=1e-12)      
        return result 
    
    def preprocess(self): 
        # repeatedly used constants 
        self.vec1 = np.ones(self.RKsteps)
        self.Ainv = np.linalg.inv(self.A)

        # not time dependent matrices Bk  
        self.Bk = np.empty(((self.m01, self.RKsteps, self.RKsteps)))
        for k in range(self.m01): 
            self.Bk[k] = np.linalg.inv(np.eye(self.RKsteps) + self.exponents01[k] * self.h * self.A)
        
        self.Cj = np.empty(((self.m12, self.RKsteps, self.RKsteps)))
        for j in range(self.m12): 
            self.Cj[j] = np.linalg.inv(np.eye(self.RKsteps) + self.exponents12[j] * self.h * self.A)
        
        # not time dependent matrix D 
        tmpD = np.zeros((self.RKsteps, self.RKsteps))
        for k in range(self.m01): 
            tmpD += self.weights01[k] * self.Bk[k]
        self.D = np.matmul(np.matmul(self.A, tmpD), self.Ainv)

        # not time dependent matrix E
        tmpE = np.zeros((self.RKsteps, self.RKsteps))
        for j in range(self.m12): 
            tmpE += self.weights12[j] * self.Cj[j] 
        self.E = np.matmul(np.matmul(self.A, tmpE), self.Ainv)
            
        # not time dependent matrix F
        self.F = np.matmul(self.E, self.Ainv) / self.h 
        
    def selectRK(self, scheme): 
        schemes = {
            "IE"  : self.setImplicitEuler, 
            "RIIA1" : self.setImplicitEuler, 
            "RIIA2" : self.setRadauIIA2, 
            "RIIA3" : self.setRadauIIA3, 
        }
        return schemes[scheme]()

    def setImplicitEuler(self): 
        self.A = np.array([[1]])
        self.b = np.array([1])
        self.c = np.array([1])
        self.RKsteps = 1 

    def setRadauIIA2(self): 
        self.A = np.array([[5/12, -1/12], 
                           [3/4, 1/4]])
        self.b = np.array([3/4, 1/4])
        self.c = np.array([1/3, 1])
        self.RKsteps = 2

    def setRadauIIA3(self): 
        sq6 = np.sqrt(6)

        self.c = np.asarray([
            (4 - sq6) / 10, 
            (4 + sq6) / 10, 
            1
        ])

        self.b = np.asarray([
            (16 - sq6) / 36, 
            (16 + sq6) / 36, 
            1. / 9.  
        ]) 

        self.A = np.zeros((3,3))
        self.A[0,0] = ( 88 -   7 * sq6) /  360
        self.A[0,1] = (296 - 169 * sq6) / 1800
        self.A[0,2] = ( -2 +   3 * sq6) /  225
        
        self.A[1,0] = (296 + 169 * sq6) / 1800
        self.A[1,1] = ( 88 +   7 * sq6) /  360
        self.A[1,2] = ( -2 -   3 * sq6) /  225
 
        self.A[2,:] = self.b

        self.RKsteps = 3 
        
# sufficiently smooth for high order convergence 

def exampleFirstMode(): 
    a = 12
    b = 15
    exponents01 = np.array([a])
    weights01 = np.array([b])
    exponents12 = np.array([0])
    weights12 = np.array([0])
    def rhs(t,y): 
        return 3 * b * (t**2/a - 2*t/a**2 - 2*np.exp(-a*t)/a**3 + 2/a**3)
    def sol(t): 
        return t**3
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "FirstMode"+str(a)+"_"+str(b) 

def exampleSecondMode(): 
    a = 1
    b = 1
    exponents01 = np.array([0])
    weights01 = np.array([0])
    exponents12 = np.array([a])
    weights12 = np.array([b])

    def rhs(t,y): 
        return 12 * b * (2 - 2*np.exp(-a*t) -2*a*t + a**2*t**2) / a**3
    def sol(t): 
        return t**4 
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "SecondMode"+str(a)+"_"+str(b) 

def exampleSingleModes(): 
    a = 1
    b = 1
    c = 1
    d = 1
    exponents01 = np.asarray([a])
    weights01 = np.asarray([b])
    exponents12 = np.asarray([c])
    weights12 = np.asarray([d])
    def rhs(t,y): 
        mode1 = 4 * b * (6 * np.exp(-a*t) -6 + 6*a*t - 3 * a**2 * t**2 + a**3 * t**3)/a**4
        mode2 = 12 * d * (2 - 2*np.exp(-c*t) -2*c*t + c**2*t**2) / c**3
        return mode1 + mode2 
    def sol(t): 
        return t**4 
    
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "SingleModes"

def exampleSingleModesFofTandU(): 
    a = 1
    b = 2
    c = 3
    d = 4
    exponents01 = np.asarray([a])
    weights01 = np.asarray([b])
    exponents12 = np.asarray([c])
    weights12 = np.asarray([d])
    def rhs(t,y): 
        mode1 = 4 * b * (6 * np.exp(-a*t) -6 + 6*a*t - 3 * a**2 * t**2 + a**3 * t**3)/a**4
        mode2 = 12 * d * (2 - 2*np.exp(-c*t) -2*c*t + c**2*t**2) / c**3
        return mode1 + mode2 - 10 * np.sin(20 * y**2) + 10 * np.sin(20 * t**8) 
    
    def sol(t): 
        return t**4 
    
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "SingleModes"

def exampleGamma6(): 
    file01 = 'data/RA_DO/distGamma6_m_69_AAAtol_1.0e-30.csv'
    data01 = np.genfromtxt(file01, delimiter=',')
    weights01 = data01[0,:]
    exponents01 = -data01[1,:] 

    file12 = 'data/RA_DO/distGamma5_m_68_AAAtol_1.0e-30.csv'
    data12 = np.genfromtxt(file12, delimiter=',')
    weights12 = data12[0,:]
    exponents12 = -data12[1,:] 

    def rhs(t,y): 
        return np.nan_to_num((t ** 5 - t**3) / np.log(t), nan=2) 
    
    def sol(t): 
        return t**5  
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "Gamma6"

def exampleGamma4Sinh(): 
    file01 = 'data/RADO/distGamma4sinh_QF_15_1.0e-50_AAA_1.0e-50_zmax_1000_nz_1001.csv'
    data01 = np.genfromtxt(file01, delimiter=',')
    weights01 = data01[0,:]
    exponents01 = -data01[1,:] 

    file12 = 'data/RADO/distGamma3sinh_QF_15_1.0e-50_AAA_1.0e-50_zmax_1000_nz_1001.csv'
    data12 = np.genfromtxt(file12, delimiter=',')
    weights12 = data12[0,:]
    exponents12 = -data12[1,:] 

    def rhs(t,y): 
        return 6 * t * (t**2 - np.cosh(2) - np.sinh(2) * np.log(t)) / (np.log(t) ** 2 - 1)
        
    def sol(t): 
        return t**3 
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "Gamma6"

def exampleExpGamma6(): 
    file01 = 'data/RA_DO/distExpGamma6_m_103_AAAtol_1.0e-45.csv'
    data01 = np.genfromtxt(file01, delimiter=',')
    weights01 = data01[0,:]
    exponents01 = -data01[1,:] 

    file12 = 'data/RA_DO/distExpGamma5_m_103_AAAtol_1.0e-45.csv'
    data12 = np.genfromtxt(file12, delimiter=',')
    weights12 = data12[0,:]
    exponents12 = -data12[1,:] 

    def rhs(t,y): 
        return 120 * (np.exp(1)**2 *t**5 - t**3) / (np.exp(1)**2 * (np.log(t) + 1))   # - 100 * np.sin(20 * (t**3-y)) + 100 * np.sin(20 * (t**3 - t**5)) 
    
    def sol(t): 
        return t**5 
    
    return 0,0, weights01, exponents01, weights12, exponents12, rhs, sol, "ExpGamma6"


def exampleGamma4(): 
    file01 = 'data/RA_DO/distGamma4_m_102_AAAtol_1.0e-45.csv'
    data01 = np.genfromtxt(file01, delimiter=',')
    weights01 = data01[0,:]
    exponents01 = -data01[1,:] 

    file12 = 'data/RA_DO/distGamma3_m_102_AAAtol_1.0e-45.csv'
    data12 = np.genfromtxt(file12, delimiter=',')
    weights12 = data12[0,:]
    exponents12 = -data12[1,:] 

    def rhs(t,y): 
        return (6*t**3 + 6*t -4) / np.log(t) + (6-10*t) / np.log(t)**2 + (4*t-4) / np.log(t)**3
    
    def sol(t): 
        return t**3 + 2*t + 4  

    u0 = 4.
    v0 = 2.
    
    return u0,v0, weights01, exponents01, weights12, exponents12, rhs, sol, "ExpGamma6"
 
def exampleSettings(name): 
    example_dict = {
        "FirstMode" : exampleFirstMode, 
        "SecondMode" : exampleSecondMode, 
        "SingleModes" : exampleSingleModes,
        "exampleSingleModesFofTandU" : exampleSingleModesFofTandU,  
        "Gamma6" : exampleGamma6, 
        "Gamma4" : exampleGamma4, 
        "Gamma4Sinh" : exampleGamma4Sinh,  
        "ExpGamma6" : exampleExpGamma6
    }     
    return example_dict[name]()

if __name__ == "__main__":
    nphyper = np.vectorize( hyper )
    nphyper.excluded.add(0)
    nphyper.excluded.add(1)

    ########################################################################
    u0, v0, weights01, exponents01, weights12, exponents12, rhs, sol, example = exampleSettings("SingleModes")

    schemes = ["RIIA1" , "RIIA2", "RIIA3"]
    schemes = ["RIIA2"]
    orders = {}
    orders["RIIA1"] = 1 
    orders["RIIA2"] = 3 
    orders["RIIA3"] = 5 

    cmap = plt.get_cmap("Dark2")

    # hs = np.geomspace(0.1,0.1,1)
    hs = np.geomspace(0.1,0.0001,7)
    ihs = 1/hs

    fig_solution, ax_solution = plt.subplots(1,3,figsize=(15,5), dpi=100)
    fig_modes, ax_modes = plt.subplots(1,3,figsize=(15,5), dpi=100)
    fig_error, ax_error = plt.subplots(1,3,figsize=(15,5), dpi=100)

    fig_L1_convergence, ax_L1_convergence = plt.subplots(1,1,figsize=(5,5), dpi=100)    
    fig_Linf_convergence, ax_Linf_convergence = plt.subplots(1,1,figsize=(5,5), dpi=100)    
    fig_LT_convergence, ax_LT_convergence = plt.subplots(1,1,figsize=(5,5), dpi=100)    

    for idx, scheme in enumerate(schemes):
        print(scheme)

        L1error = []
        Linferror = []
        LTerror = []

        for h in hs: 
            print(h)
            COFDE = FDE02RK(weights01=weights01, exponents01=exponents01, 
                            weights12=weights12, exponents12=exponents12, 
                            f=rhs, t0=0, u0=np.array([u0]), v0=np.array([v0]), tf=1.0, h=h, scheme=scheme)
            times, numsol = COFDE.solve()
            error = np.abs(numsol[:,0] - sol(times))
            
            ax_solution[idx].plot(times,numsol[:,0], label='{:0.3e}'.format(h))
            # ax_solution[idx].plot(times,numsol[:,1], label='{:0.3e}'.format(h))
            ax_modes[idx].plot(times,np.sum(numsol[:,2:], axis=1), label='{:0.3e}'.format(h))
            ax_error[idx].plot(times,error, label='{:0.3e}'.format(h))

            L1error.append(np.mean(error))
            Linferror.append(np.max(error))
            LTerror.append(error[-1])

        fig_solution.suptitle(example + " solution")
        ax_solution[idx].plot(times,sol(times), label='solution', color="black", linewidth=2)
        # ax_solution[idx].plot(times,times**2, label='derivative', color="black", linewidth=2)
        ax_solution[idx].set_title(scheme)
        ax_solution[idx].legend()

        fig_modes.suptitle(example + " sum over all modes")
        ax_modes[idx].plot(times,rhs(times,sol(times)), label='rhs', color="black", linewidth=2)
        ax_modes[idx].set_title(scheme)
        ax_modes[idx].legend()

        fig_error.suptitle(example + " absolute error")
        ax_error[idx].set_title(scheme)
        ax_error[idx].legend()
        ax_error[idx].set_yscale("log")
        ax_error[idx].set_ylim([1e-25,1e2])

        ax_L1_convergence.set_title(example + " L1 error")
        ax_L1_convergence.plot(ihs, L1error, label=scheme, marker="s", color=cmap(idx/5-1e-12))
        ax_L1_convergence.plot(ihs, (hs/hs[0])**orders[scheme] * L1error[0], linestyle="--", color=cmap(idx/5-1e-12))
        ax_L1_convergence.plot(ihs, (hs/hs[0])**(orders[scheme]-1) * L1error[0], linestyle=":", color=cmap(idx/5-1e-12))
        ax_L1_convergence.set_xscale("log")
        ax_L1_convergence.set_yscale("log")
        ax_L1_convergence.set_ylim([1e-20,1e2])
        ax_L1_convergence.legend()
        ax_L1_convergence.grid('on')
        
        ax_Linf_convergence.set_title(example + " Linf error")
        ax_Linf_convergence.plot(ihs, Linferror, label=scheme, marker="s", color=cmap(idx/5-1e-12))
        ax_Linf_convergence.plot(ihs, (hs/hs[0])**orders[scheme] * Linferror[0], linestyle="--", color=cmap(idx/5-1e-12))
        ax_Linf_convergence.plot(ihs, (hs/hs[0])**(orders[scheme]-1) * Linferror[0], linestyle=":", color=cmap(idx/5-1e-12))
        ax_Linf_convergence.set_xscale("log")
        ax_Linf_convergence.set_yscale("log")
        ax_Linf_convergence.set_ylim([1e-20,1e2])
        ax_Linf_convergence.legend()
        ax_Linf_convergence.grid('on')

        ax_LT_convergence.set_title(example + " LT error")
        ax_LT_convergence.plot(ihs, LTerror, label=scheme, marker="s", color=cmap(idx/5-1e-12))
        ax_LT_convergence.plot(ihs, (hs/hs[0])**orders[scheme] * LTerror[0], linestyle="--", color=cmap(idx/5-1e-12))
        ax_LT_convergence.plot(ihs, (hs/hs[0])**(orders[scheme]-1) * LTerror[0], linestyle=":", color=cmap(idx/5-1e-12))
        ax_LT_convergence.set_xscale("log")
        ax_LT_convergence.set_yscale("log")
        ax_LT_convergence.set_ylim([1e-20,1e2])
        ax_LT_convergence.legend()
        ax_LT_convergence.grid('on')
    
    # fig_solution.savefig("Collocation_methods_" + example + "_solution.pdf", bbox_inches="tight")
    # fig_modes.savefig("Collocation_methods_" + example + "_modes.pdf", bbox_inches="tight")
    # fig_error.savefig("Collocation_methods_" + example + "_error.pdf", bbox_inches="tight")

    # fig_L1_convergence.savefig("Collocation_methods_" + example + "_L1convergence.pdf", bbox_inches="tight")
    # fig_Linf_convergence.savefig("Collocation_methods_" + example + "_Linfconvergence.pdf", bbox_inches="tight")
    # fig_LT_convergence.savefig("Collocation_methods_" + example + "_LTconvergence.pdf", bbox_inches="tight")
    plt.show()

