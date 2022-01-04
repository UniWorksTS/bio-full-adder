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
    
   
def not_not_yes_or(A, B, C, N_A, N_B, N_C, L_A, L_B, out, r_X, params_not, params_yes):
    # NOT A or NOT B or C

    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)

    state_not_B = L_B, out, B, N_B
    dL_B, dd = not_cell_wrapper(state_not_B, params_not)
    d_out += dd
    dN_B = population(N_B, r_X)

    state_yes_C = out, C, N_C
    d_out += yes_cell_wrapper(state_yes_C, params_yes)
    dN_C = population(N_C, r_X)


    return dN_A, dN_B, dN_C, dL_A, dL_B, d_out
    
    
def not_yes_yes_or(A, B, C, N_A, N_B, N_C, L_A, out, r_X, params_not, params_yes):
    # NOT A or NOT B or C

    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)

    state_yes_B = out, B, N_B
    d_out += yes_cell_wrapper(state_yes_B, params_yes)
    dN_B = population(N_B, r_X)

    state_yes_C = out, C, N_C
    d_out += yes_cell_wrapper(state_yes_C, params_yes)
    dN_C = population(N_C, r_X)


    return dN_A, dN_B, dN_C, dL_A, d_out
    
    
def nand3(A, B, C, N_A, N_B, N_C, L_A, L_B, out, r_X, params_not, params_yes):

    d_out = 0
    
    
    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)

    state_not_B = L_B, out, B, N_B
    dL_B, dd = not_cell_wrapper(state_not_B, params_not)
    d_out += dd
    dN_B = population(N_B, r_X)
    
    state_not_C = L_C, out, C, N_C
    dL_C, dd = not_cell_wrapper(state_not_C, params_not)
    d_out += dd
    dN_C = population(N_C, r_X)

    return dN_A, dN_B, dN_C, dL_A, dL_B, dL_C, d_out
    
    
def not_or4(A, B, C, D, N_A, N_B, N_C, N_D, L_A, L_B, L_C, L_D, out, r_X, params_not):
    # NOT A or NOT B or NOT C or NOT D

    d_out = 0

    state_not_A = L_A, out, A, N_A
    dL_A, dd = not_cell_wrapper(state_not_A, params_not)
    d_out += dd
    dN_A = population(N_A, r_X)

    state_not_B = L_B, out, B, N_B
    dL_B, dd = not_cell_wrapper(state_not_B, params_not)
    d_out += dd
    dN_B = population(N_B, r_X)

    state_not_C = L_C, out, C, N_C
    dL_C, dd = not_cell_wrapper(state_not_C, params_not)
    d_out += dd
    dN_C = population(N_C, r_X)

    state_not_D = L_D, out, D, N_D
    dL_D, dd = not_cell_wrapper(state_not_D, params_not)
    d_out += dd
    dN_D = population(N_D, r_X)

    return dN_A, dN_B, dN_C, dN_D, dL_A, dL_B, dL_C, dL_D, d_out
    
    
    
def or4(A, B, C, D, N_A, N_B, N_C, N_D, out, r_X, params_yes):
    # YES x OR YES y = OR

    d_out = 0

    state_yes_A = out, A, N_A
    d_out += yes_cell_wrapper(state_yes_A, params_yes)
    dN_A = population(N_A, r_X)

    state_yes_B = out, B, N_B
    d_out += yes_cell_wrapper(state_yes_B, params_yes)
    dN_B = population(N_B, r_X)
    
    state_yes_C = out, C, N_C
    d_out += yes_cell_wrapper(state_yes_C, params_yes)
    dN_C = population(N_C, r_X)
    
    state_yes_D = out, D, N_D
    d_out += yes_cell_wrapper(state_yes_D, params_yes)
    dN_D = population(N_D, r_X)

    return dN_A, dN_B, dN_C, dN_D, d_out
    
    
def sum_model(state, params):
    L_A1, L_A2, L_A3, L_A4, L_B1, L_B2, L_B3, L_B4, L_Ci1, L_Ci2, L_Ci3, L_Ci4, L_AnotBnotC, L_notABnotC, L_notAnotBC, L_ABnotC, a, b, ci, sumout, AnotBnotC, notABnotC, notAnotBC, ABnotC = state

    stateA1 = L_A1, AnotBnotC, a, 1, 1
    stateB1 = L_B1, AnotBnotC, b, 1, 1
    stateC1 = L_Ci1, AnotBnotC, ci, 1, 1
    stateAnotBnotC = L_AnotBnotC, sumout, AnotBnotC, 1, 1

    stateA2 = L_A2, notABnotC, a, 1, 1
    stateB2 = L_B2, notABnotC, b, 1, 1
    stateC2 = L_Ci2, notABnotC, ci, 1, 1
    statenotABnotC = L_notABnotC, sumout, notABnotC, 1, 1

    stateA3 = L_A3, notAnotBC, a, 1, 1
    stateB3 = L_B3, notAnotBC, b, 1, 1
    stateC3 = L_Ci3, notAnotBC, ci, 1, 1
    statenotAnotBC = L_notAnotBC, sumout, notAnotBC, 1, 1

    stateA4 = L_A4, ABnotC, a, 1, 1
    stateB4 = L_B4, ABnotC, b, 1, 1
    stateC4 = L_Ci4, ABnotC, ci, 1, 1
    stateABnotC = L_ABnotC, sumout, ABnotC, 1, 1

    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a

    dL_A1_dt, dAnotBnotC1_dt = not_cell(stateA1, params_A)
    dL_B1_dt, dAnotBnotC2_dt = not_cell(stateB1, params_A)
    dL_C1_dt, dAnotBnotC3_dt = not_cell(stateC1, params_A)
    dAnotBnotC_dt = dAnotBnotC1_dt + dAnotBnotC2_dt + dAnotBnotC3_dt
    dl_AnotBnotC_dt, dsumout1_dt = not_cell(stateAnotBnotC, params_A)

    dL_A2_dt, dnotABnotC1_dt = not_cell(stateA2, params_A)
    dL_B2_dt, dnotABnotC2_dt = not_cell(stateB2, params_A)
    dL_C2_dt, dnotABnotC3_dt = not_cell(stateC2, params_A)
    dnotABnotC_dt = dnotABnotC1_dt + dnotABnotC2_dt + dnotABnotC3_dt
    dl_notABnotC_dt, dsumout2_dt = not_cell(statenotABnotC, params_A)

    dL_A3_dt, dnotAnotBC1_dt = not_cell(stateA3, params_A)
    dL_B3_dt, dnotAnotBC2_dt = not_cell(stateB3, params_A)
    dL_C3_dt, dnotAnotBC3_dt = not_cell(stateC3, params_A)
    dnotAnotBC_dt = dnotAnotBC1_dt + dnotAnotBC2_dt + dnotAnotBC3_dt
    dl_notAnotBC_dt, dsumout3_dt = not_cell(statenotAnotBC, params_A)

    dL_A4_dt, dABnotC1_dt = not_cell(stateA4, params_A)
    dL_B4_dt, dABnotC2_dt = not_cell(stateB4, params_A)
    dL_C4_dt, dABnotC3_dt = not_cell(stateC4, params_A)
    dABnotC_dt = dABnotC1_dt + dABnotC2_dt + dABnotC3_dt
    dl_ABnotC_dt, dsumout4_dt = not_cell(stateABnotC, params_A)



    dsumout_dt = dsumout1_dt + dsumout2_dt + dsumout3_dt + dsumout4_dt
        
    return np.array([dL_A1_dt, dL_A2_dt, dL_A3_dt, dL_A4_dt, dL_B1_dt, dL_B2_dt, dL_B3_dt, dL_B4_dt, dL_C1_dt, dL_C2_dt, dL_C3_dt, dL_C4_dt, dl_AnotBnotC_dt, dl_notABnotC_dt, dl_ABnotC_dt, dl_notAnotBC_dt, 0, 0, 0, dsumout_dt, dAnotBnotC_dt, dnotABnotC_dt, dnotAnotBC_dt, dABnotC_dt])

def sum_model_ODE(T, state, params):
    return sum_model(state, params)

def carryOut_model(state, params):
    L_A1, L_A2, L_B1, L_B2, L_Ci1, L_Ci2, L_AB, L_AC, L_BC, a, b, ci, cout, ab, ac, bc = state

    stateA1 = L_A1, ab, a, 1, 1
    stateB1 = L_B1, ab, b, 1, 1
    stateAB = L_AB, cout, ab, 1, 1
    
    stateA2 = L_A2, ac, a, 1, 1
    stateC1 = L_Ci1, ac, ci, 1, 1
    stateAC = L_AC, cout, ac, 1, 1

    stateB2 = L_B2, bc, b, 1, 1
    stateC2 = L_Ci2, bc, ci, 1, 1
    stateBC = L_BC, cout, bc, 1, 1

    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a

    dL_A1_dt, dab1_dt = not_cell(stateA1, params_A)
    dL_B1_dt, dab2_dt = not_cell(stateB1, params_A)
    dab_dt = dab1_dt + dab2_dt
    dl_AB_dt, dcout1_dt = not_cell(stateAB, params_A)

    dL_A2_dt, dac1_dt = not_cell(stateA2, params_A)
    dL_C1_dt, dac2_dt = not_cell(stateC1, params_A)
    dac_dt = dac1_dt + dac2_dt
    dl_AC_dt, dcout2_dt = not_cell(stateAC, params_A)

    dL_B2_dt, dbc1_dt = not_cell(stateB2, params_A)
    dL_C2_dt, dbc2_dt = not_cell(stateC2, params_A)
    dbc_dt = dbc1_dt + dbc2_dt
    dl_BC_dt, dcout3_dt = not_cell(stateBC, params_A)
    
    dcout_dt = dcout1_dt + dcout2_dt + dcout3_dt
    #print("out : " + str(dcout_dt))
        
    return np.array([dL_A1_dt, dL_A2_dt, dL_B1_dt, dL_B2_dt, dL_C1_dt, dL_C2_dt, dl_AB_dt, dl_AC_dt, dl_BC_dt, 0, 0, 0, dcout_dt, dab_dt, dac_dt, dbc_dt])

def carryOut_model_ODE(T, state, params):
    return carryOut_model(state, params)
    
params = delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, rho_x #delta_y, rho_x, rho_y, r_X 

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
#print(not_model(Y0, T[0], params))

#print(nand(1,0,1,1,1,1, T[0],r_X, params))
