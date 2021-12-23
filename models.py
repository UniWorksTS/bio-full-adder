## models of cel
from parameters import *
import numpy as np

def yes_cell(state, params):
    x, y, N_X, N_Y = state
    gamma_x, n_y, theta_x, delta_x, rho_x = params

    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    dx_dt = N_X * gamma_x * (y ** n_y)/(1 + (theta_x*y)**n_y ) - N_Y * (delta_x * x) - rho_x * x
    
    return dx_dt
    
    
    
def not_cell(state, params):
    L_X, x, y, N_X, N_Y = state
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x = params


    # presume that the molecules are degraded in the same strain as they are produced
    N_Y = N_X


    f = gamma_L_X * (y ** n_y)/(1 + (theta_L_X*y)**n_y )
    dL_X_dt = N_X * (f - delta_L * L_X)

    dx_dt = N_X * (eta_x * (1/(1+ (omega_x*L_X)**m_x))) - N_Y * (delta_x * x) - rho_x * x

    return dL_X_dt, dx_dt
    
    
def population(state, params):
    N = state
    r = params

    dN = r * N * (1 - N)    

    return dN
    
    
# L_A ... intermediate
# a ... out
# b ... in
# N_A ... number of cells
def not_cell_wrapper(state, params):
    L_A, a, b, N_A = state

    state_A = L_A, a, b, N_A, N_A
    params_A = params

    return not_cell(state_A, params_A)


# a ... out
# b ... in
# N_A ... number of cells
def yes_cell_wrapper(state, params):
    a, b, N_A = state

    state_A = a, b, N_A, N_A
    params_A = params

    return yes_cell(state_A, params_A)
    
    
def not_not_yes_or(A, B, C, N_A, N_B, N_C, L_A, L_B, out, r_X, params_not, params_yes):
    # NOT A or NOT B or C

    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)
    print(d_out)

    state_not_B = L_B, out, B, N_B
    dL_B, dd = not_cell_wrapper(state_not_B, params_not)
    d_out += dd
    dN_B = population(N_B, r_X)
    print(d_out)
    state_yes_C = out, C, N_C
    d_out += yes_cell_wrapper(state_yes_C, params_yes)
    dN_C = population(N_C, r_X)

    return dN_A, dN_B, dN_C, dL_A, dL_B, d_out
    
    
def not_yes_yes_or(A, B, C, N_A, N_B, N_C, L_A, out, r_X, params_not, params_yes):
    # NOT A or B or C

    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)
    print(d_out)

    state_yes_B = out, B, N_B
    d_out += yes_cell_wrapper(state_yes_B, params_yes)
    dN_B = population(N_B, r_X)
    print(d_out)
    
    state_yes_C = out, C, N_C
    d_out += yes_cell_wrapper(state_yes_C, params_yes)
    dN_C = population(N_C, r_X)

    return dN_A, dN_B, dN_C, dL_A, d_out
    
def or2(A, B, N_A, N_B, out, r_X, params_yes):
    # YES x OR YES y = OR

    d_out = 0

    state_yes_A = out, A, N_A
    d_out += yes_cell_wrapper(state_yes_A, params_yes)
    dN_A = population(N_A, r_X)

    state_yes_B = out, B, N_B
    d_out += yes_cell_wrapper(state_yes_B, params_yes)
    dN_B = population(N_B, r_X)

    return dN_A, dN_B, d_out


def nand(A, B, N_A, N_B, L_A, L_B, out, r_X, params_not):
    # NOT A OR NOT B = NAND
    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)

    state_not_B = L_B, out, B, N_B
    dL_B, dd = not_cell_wrapper(state_not_B, params_not)
    d_out += dd
    dN_B = population(N_B, r_X)

    return dN_A, dN_B, dL_A, dL_B, d_out
   
    
    
def not_model(state, T, params):
    L_A, a, b, N_A = state

    delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, delta_b, rho_a, rho_b, r_A = params

    state_not = L_A, a, b, N_A
    params_not = delta_L, gamma_L_A, n_b, theta_L_A, eta_a, omega_a, m_a, delta_a, rho_a
    dL_A_dt, da_dt = not_cell_wrapper(state_not, params_not)
    
    db_dt = 0#- N_A * delta_b * b - rho_b * b

    dN_A_dt = population(N_A, r_A)

    return np.array([dL_A_dt, da_dt, db_dt, dN_A_dt])

def yes_model(state, T, params):
    a, b, N_A = state
    
    gamma_a, n_b, theta_a, delta_a, delta_b, rho_a, rho_b, r_A = params

    state_yes = a, b, N_A
    params_yes = gamma_a, n_b, theta_a, delta_a, rho_a
    da_dt = yes_cell_wrapper(state_yes, params_yes)
    
    db_dt = 0 #- N_A * delta_b * b - rho_b * b

    dN_A_dt = population(N_A, r_A)

    return np.array([da_dt, db_dt, dN_A_dt])
    
    
def sum_model(state, T, params):

    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x, gamma_x, theta_x, r_X = params
    params_yes = gamma_x, n_y, theta_x, delta_x, rho_x
    params_not = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x
    out = 0.5
    N_A, N_B, N_C, L_A, L_B = 1,1,1,1,1
    #
    A, B, C_IN = state[:3]
    
    
    m_4_3 = nand(A,B, N_A, N_B, L_A, L_B, out, r_X, params_not)
    print(m_4_3)

    
    
    
params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, r_X 


# simulation parameters
t_end = 1500
N = t_end

# Y = L_A, a, b, N_A
Y0 = np.zeros(4)
Y0[1] = 1
Y0[2] = 0 # b
Y0[3] = 1 # N_A

T = np.linspace(0, t_end, N)



t1 = t_end
dt = t_end/N
T = np.arange(0,t1+dt,dt)
Y = np.zeros([1+N,4])
Y[0,:] = Y0 
print(not_model(Y0, T[0], params))