from scipy.integrate import ode
import matplotlib.pyplot as plt

from models import *
from parameters import *

params = [delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B]

params[-4] = 0
params[-3] = 0

# simulation parameters
t_end = 100
N = t_end

for a in range(2):
    for b in range(2):
        for c in range(2):
            fig, axs = plt.subplots(2)
            fig.suptitle('simulation results')
            Y0 = np.zeros(16)
            Y0[-7] = a * 10#A
            Y0[-6] = b * 10#B
            Y0[-5] = c * 10#Ci
            #Y0[-4] = 1 #Co
            T = np.linspace(0, t_end, N)

            t1 = t_end
            dt = t_end/N
            T = np.arange(0,t1+dt,dt)
            Y = np.zeros([1+N,16])
            Y[0,:] = Y0

            r = ode(carryOut_model_ODE).set_integrator('zvode', method='bdf')
            r.set_initial_value(Y0, T[0]).set_f_params(params)

            i = 1
            while r.successful() and r.t < t1:
                Y[i,:] = r.integrate(r.t+dt)
                i += 1
            
            AC = Y[:,-7]
            BC = Y[:,-6]
            CinC = Y[:,-5]
            CoutC = Y[:,-4]
            ABC = Y[:,-3]
            ACC = Y[:,-2]
            BCC = Y[:,-1]
            axs[0].plot(T, AC)
            axs[0].plot(T, BC)
            axs[0].plot(T, CinC)
            axs[0].plot(T, CoutC)
            #axs[0].plot(T, ABC)
            #axs[0].plot(T, ACC)
            #axs[0].plot(T, BCC)
            axs[0].legend(["A", "B", "Cin", "Cout"])

            Z0 = np.zeros(18)
            Z0[-8] = a * 10#A
            Z0[-7] = b * 10#B
            Z0[-6] = c * 10#Ci
            #Y0[-5] = 1 #sum
            T = np.linspace(0, t_end, N)

            t1 = t_end
            dt = t_end/N
            T = np.arange(0,t1+dt,dt)
            Z = np.zeros([1+N,18])
            Z[0,:] = Z0

            q = ode(sum_model_ODE).set_integrator('zvode', method='bdf')
            q.set_initial_value(Z0, T[0]).set_f_params(params)

            i = 1
            while q.successful() and q.t < t1:
                Z[i,:] = q.integrate(q.t+dt)
                i += 1
            
            AS = Z[:,-7]
            BS = Z[:,-6]
            CinS = Z[:,-5]
            CoutS = Z[:,-4]
            ABS = Z[:,-3]
            ACS = Z[:,-2]
            BCS = Z[:,-1]
            axs[1].plot(T, AS)
            axs[1].plot(T, BS)
            axs[1].plot(T, CinS)
            axs[1].plot(T, CoutS)
            #axs[1].plot(T, ABS)
            #axs[1].plot(T, ACS)
            #axs[1].plot(T, BCS)
            axs[1].legend(["A", "B", "Cin", "sum"])
            plt.show()