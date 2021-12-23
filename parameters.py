# see Supplementary Table 1 in Urrios 2016

gamma_A = 0.615 #nM/min
gamma_B = 0.495 #nM/min 
mu_A = 2 #
mu_B = 2 # 
omega_a = 1550 #nM-1
omega_b = 1550 #nM-1
eta_a = 0.0369#0.0369 #nM/min 
eta_b = 0.162#0.162 #nM/min 
r_A = 0.07 #10
r_B = 0.07 #10
m_a = 2 
m_b = 2
delta_L = 0.15 #min-1
delta_a = 0.05 #min-1
delta_b = 0.023 #min-1
theta_A = 0.26 #nM-1
theta_B = 0.167 #nM-1
n_a = 0.9
n_b = 1.2
rho_a = 5#5 #min-1
rho_b = 5#5 #min-1

gamma_a = gamma_A
gamma_b = gamma_B
theta_a = theta_A
theta_b = theta_B


scenario = "b"

if scenario == "a":
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x = delta_L, gamma_A, n_a, theta_A, eta_a, omega_a, m_a, delta_a, delta_a, rho_a, rho_a, gamma_a, theta_a
    r_X, r_Y = r_A, r_A
elif scenario == "b":
    delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x = delta_L, gamma_B, n_b, theta_B, eta_b, omega_b, m_b, delta_b, delta_b, rho_b, rho_b, gamma_b, theta_b
    r_X, r_Y = r_B, r_B


print(delta_L, gamma_L_X, n_y, theta_L_X, eta_x, omega_x, m_x, delta_x, delta_y, rho_x, rho_y, gamma_x, theta_x, r_X, r_Y)