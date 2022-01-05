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
    
    
    
def sum_model(state, params):
    L_A1, L_A4, L_B2, L_B4, L_Ci3, L_Ci4, L_AnotBnotC, L_notABnotC, L_notAnotBC, L_ABC, a, b, ci, sumout, AnotBnotC, notABnotC, notAnotBC, ABC = state

    
    stateA1 = L_A1, AnotBnotC, a, 1, 1
    stateB1 = AnotBnotC, b, 1, 1
    stateC1 = AnotBnotC, ci, 1, 1
    stateAnotBnotC = L_AnotBnotC, sumout, AnotBnotC, 1, 1

    stateA2 = notABnotC, a, 1, 1
    stateB2 = L_B2, notABnotC, b, 1, 1
    stateC2 = notABnotC, ci, 1, 1
    statenotABnotC = L_notABnotC, sumout, notABnotC, 1, 1

    stateA3 = notAnotBC, a, 1, 1
    stateB3 = notAnotBC, b, 1, 1
    stateC3 = L_Ci3, notAnotBC, ci, 1, 1
    statenotAnotBC = L_notAnotBC, sumout, notAnotBC, 1, 1

    stateA4 = L_A4, ABC, a, 1, 1
    stateB4 = L_B4, ABC, b, 1, 1
    stateC4 = L_Ci4, ABC, ci, 1, 1
    stateABC = L_ABC, sumout, ABC, 1, 1

    delta_L, gamma_A, gamma_B, n_a, n_b, theta_A, theta_B, eta_a, eta_b, omega_a, omega_b, m_a, m_b, delta_a, delta_b, rho_a, rho_b, r_A, r_B = params

    params_A = delta_L, gamma_A, n_b, theta_A, eta_a, omega_a, m_a, delta_a, rho_a
    params_B = gamma_B, n_a, theta_B, delta_b, rho_b


    dL_A1_dt, dAnotBnotC1_dt = not_cell(stateA1, params_A)
    dAnotBnotC2_dt = yes_cell(stateB1, params_B)
    dAnotBnotC3_dt = yes_cell(stateC1, params_B)
    dAnotBnotC_dt = dAnotBnotC1_dt + dAnotBnotC2_dt + dAnotBnotC3_dt
    dl_AnotBnotC_dt, dsumout1_dt = not_cell(stateAnotBnotC, params_A)

    dnotABnotC1_dt = yes_cell(stateA2, params_B)
    dL_B2_dt, dnotABnotC2_dt = not_cell(stateB2, params_A)
    dnotABnotC3_dt = yes_cell(stateC2, params_B)
    dnotABnotC_dt = dnotABnotC1_dt + dnotABnotC2_dt + dnotABnotC3_dt
    dl_notABnotC_dt, dsumout2_dt = not_cell(statenotABnotC, params_A)

    dnotAnotBC1_dt = yes_cell(stateA3, params_B)
    dnotAnotBC2_dt = yes_cell(stateB3, params_B)
    dL_C3_dt, dnotAnotBC3_dt = not_cell(stateC3, params_A)
    dnotAnotBC_dt = dnotAnotBC1_dt + dnotAnotBC2_dt + dnotAnotBC3_dt
    dl_notAnotBC_dt, dsumout3_dt = not_cell(statenotAnotBC, params_A)

    dL_A4_dt, dABC1_dt = not_cell(stateA4, params_A)
    dL_B4_dt, dABC2_dt = not_cell(stateB4, params_A)
    dL_C4_dt, dABC3_dt = not_cell(stateC4, params_A)
    dABC_dt = dABC1_dt + dABC2_dt + dABC3_dt
    dl_ABC_dt, dsumout4_dt = not_cell(stateABC, params_A)



    dsumout_dt = dsumout1_dt + dsumout2_dt + dsumout3_dt + dsumout4_dt
        
    return np.array([dL_A1_dt, dL_A4_dt, dL_B2_dt, dL_B4_dt, dL_C3_dt, dL_C4_dt, dl_AnotBnotC_dt, dl_notABnotC_dt, dl_notAnotBC_dt, dl_ABC_dt, 0, 0, 0, dsumout_dt, dAnotBnotC_dt, dnotABnotC_dt, dnotAnotBC_dt, dABC_dt])

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
