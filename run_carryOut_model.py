from scipy.integrate import ode
import matplotlib.pyplot as plt

from models import *
from parameters import *

rho_x = 0
rho_y = 0

#params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, r_X
params = delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B

# simulation parameters
t_end = 1500
N = t_end

# Y = L_A, a, b, N_A
Y0 = np.zeros(20)
Y0[-7] = 1
Y0[-6] = 1
Y0[-5] = 1
Y0[-4] = 1
Y0[-3] = 1
Y0[-2] = 1
Y0[-1] = 1
Y0[2] = 0
T = np.linspace(0, t_end, N)

t1 = t_end
dt = t_end/N
T = np.arange(0,t1+dt,dt)
Y = np.zeros([1+N,20])
Y[0,:] = Y0

r = ode(carryOut_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T[0]).set_f_params(params)

i = 1
while r.successful() and r.t < t1:
    Y[i,:] = r.integrate(r.t+dt)    # not work
    i += 1

"""
# Y = L_A, a, b, N_A
L_A = Y[:,0]
a = Y[:,1]
b = Y[:,2]
N_A = Y[:,3]

ax1 = plt.subplot(311)
ax1.plot(T,L_A)
ax1.legend(["L_A"])

ax2 = plt.subplot(312)
ax2.plot(T,a)
ax2.plot(T,b)
ax2.legend(["a", "b"])

ax3 = plt.subplot(313)
ax3.plot(T,N_A)
ax3.legend(["N_A"])

plt.show()
"""