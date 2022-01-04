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

Y0 = np.zeros(24)
Y0[-8] = 1 #A
Y0[-7] = 0 #B
Y0[-6] = 0 #Ci
#Y0[-5] = 1 #sum
T = np.linspace(0, t_end, N)

t1 = t_end
dt = t_end/N
T = np.arange(0,t1+dt,dt)
Y = np.zeros([1+N,24])
Y[0,:] = Y0

r = ode(sum_model_ODE).set_integrator('zvode', method='bdf')
r.set_initial_value(Y0, T[0]).set_f_params(params)

i = 1
while r.successful() and r.t < t1:
    Y[i,:] = r.integrate(r.t+dt)
    i += 1

A = Y[:,-8]
B = Y[:,-7]
Cin = Y[:,-6]
sum = Y[:,-5]
AnotBnotC = Y[:,-4]
notABnotC = Y[:,-3]
notAnotBC = Y[:,-2]
ABnotC = Y[:,-1]
plt.plot(T, A)
plt.plot(T, B)
plt.plot(T, Cin)
plt.plot(T, sum)
plt.plot(T, AnotBnotC)
plt.plot(T, notABnotC)
plt.plot(T, notAnotBC)
plt.plot(T, ABnotC)
plt.legend(["A", "B", "Cin", "sum", "AnotBnotC", "notABnotC","notAnotBC", "ABnotC"])
plt.show()