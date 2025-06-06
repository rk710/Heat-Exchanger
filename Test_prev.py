import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution, root_scalar
from constants import *
from NTU_method import NTU_method

#replaced N with N_tubes_per_pass
def Re_tube_calc(m_dot_2, N_tubes_per_pass):
    rho_w = constants["rho_w"]
    d_i = constants["d_i"]
    mu = constants["mu"]
#m_dot_tube edited
    m_dot_tube = m_dot_2/N_tubes_per_pass
    v_tube = m_dot_tube*4/(rho_w*np.pi*d_i**2)
    Re_tube = rho_w*v_tube*d_i/mu
    return Re_tube

def Re_sh_calc(m_dot_1, N_B, Y, N_shell_passes):
    L = constants["L"]
    d_sh = constants["d_sh"]
    d_o = constants["d_o"]
    rho_w = constants["rho_w"]
    mu = constants["mu"]
    B = L/(N_B+1)
    A_sh = d_sh/Y * (Y-d_o) * B * (1/N_shell_passes)
    v_sh = m_dot_1 / (rho_w*A_sh)
    d_chic_sh = d_sh * A_sh*4/(np.pi*d_sh**2) #Varies with number of shell passes
    Re_sh = rho_w * v_sh * d_chic_sh / mu
    return Re_sh

def heat_exchanger_pressure_drop(m_dot_1, m_dot_2, N_tubes_per_pass, N_tube_passes, N_shell_passes, N_B, Y, N_rows):
    #Hot side analysis
    rho_w = constants["rho_w"]
    d_i = constants["d_i"]
    d_noz = constants["d_noz"]
    d_sh = constants["d_sh"]
    d_o = constants["d_o"]
    L = constants ["L"]
    B = L/(N_B+1) #baffle spacing
    m_dot_tube = m_dot_2/N_tubes_per_pass
    v_tube = m_dot_tube*4/(rho_w*np.pi*d_i**2)
    Re_tube = Re_tube_calc(m_dot_2, N_tubes_per_pass)
    v_noz_2 = m_dot_2*4/(rho_w*np.pi*d_noz**2)
    f=(1.82*math.log(Re_tube, 10) - 1.64)**-2
    delta_p_tube = 0.5 * rho_w * v_tube**2 * f * N_tube_passes * L/d_i
    sigma = N_tubes_per_pass*N_tube_passes*d_i**2/d_sh**2
    K_c = 0.4*(1-0.75*sigma) + 0.1*10000/Re_tube
    K_e = (1-2*sigma+1*sigma**2) - sigma*0.1*10000/Re_tube #Obtain K_c and K_e from figure 8 in handout, function of sigma hence are just constant. Use turbulent values.
    delta_p_ends = 0.5*rho_w*v_tube**2 * (K_c + K_e) * N_tube_passes
    delta_p_noz_2 = rho_w * v_noz_2**2
    delta_p_2 = (delta_p_tube + delta_p_ends + delta_p_noz_2)/100000 #gives pressure in bar

    #Cold side analysis
    A_sh = d_sh/Y * (Y-d_o) * B * (1/N_shell_passes)
    v_sh = m_dot_1 / (rho_w*A_sh)
    Re_sh = Re_sh_calc(m_dot_1, N_B, Y, N_shell_passes)                  
    alpha = 0.2 #0.2 triangular pitch, 0.34 square
    delta_p_sh = 4*alpha*Re_sh**-0.15 *N_rows*N_shell_passes*rho_w*v_sh**2*(N_B+1) #FIX NEEDED CURRRENTLY ASSUMING 3 TUBE ROWS MAY NOT BE OPTIMAL
    v_noz_1 = m_dot_1*4/ (rho_w*np.pi*d_noz**2)
    delta_p_noz_1 = rho_w * v_noz_1**2
    R_L = 0.2655
    S_w = 0.33 * np.pi * d_sh**2 / 4
    N_CW = 0.8 * 4 #set to 1 tube in window atm, can change
    delta_p_window = N_B * R_L * (2+0.6*N_CW) * m_dot_1**2 / (2 * rho_w * A_sh * S_w)
    delta_p_1 = (delta_p_sh + delta_p_noz_1 + delta_p_window)/100000
    return delta_p_1, delta_p_2 #Return (delta_p_1, delta_p_2) in bar

def find_flow_rates(m_dot_1, m_dot_2, N_tubes_per_pass, N_tube_passes, N_shell_passes, N_B, Y, N_rows, tol=1e-3, max_iter=1000):
    rho_w = constants["rho_w"]

    for i in range(max_iter):
        delta_p_1, delta_p_2 = heat_exchanger_pressure_drop(m_dot_1, m_dot_2, N_tubes_per_pass, N_tube_passes, N_shell_passes, N_B, Y, N_rows)
        l_dot_1 = m_dot_1/(rho_w/1000)
        l_dot_2 = m_dot_2/(rho_w/1000)
        pressure_cold = float(cold_interp(l_dot_1))
        pressure_hot = float(hot_interp(l_dot_2))

        err_cold = abs(delta_p_1 - pressure_cold)
        err_hot = abs(delta_p_2 - pressure_hot)

        if err_cold < tol and err_hot < tol:
#            print(f"Converged in {i} iterations")
            break

        if pressure_cold < delta_p_1:
            m_dot_1 -= 0.001
        else:
            m_dot_1 += 0.001

        if pressure_hot < delta_p_2:
            m_dot_2 -= 0.001
        else:
            m_dot_2 += 0.001
    return m_dot_1, m_dot_2, delta_p_1, delta_p_2, pressure_cold, pressure_hot

#LMTD Function
def delta_T_lm_counter(T_cold_in, T_cold_out, T_hot_in, T_hot_out):
    dT1 = T_hot_in - T_cold_out
    dT2 = T_hot_out - T_cold_in
    #Avoid division by zero or log of zero
    if dT1 == dT2:
        return dT1
    if dT1 <= 0 or dT2 <= 0:
        return 1e-6  #Small positive value to penalize invalid regions
    return (dT1 - dT2) / np.log(dT1 / dT2)

def delta_T_lm_parallel(T_cold_in, T_cold_out, T_hot_in, T_hot_out):
    dT1 = T_hot_in - T_cold_in
    dT2 = T_hot_out - T_cold_out
    #Avoid division by zero or log of zero
    if dT1 == dT2:
        return dT1
    if dT2 <= 0:
        return 1e-6  #Small positive value to penalize invalid regions
    return (dT1 - dT2) / np.log(dT1 / dT2)

def F_1_2N(T_cold_in, T_cold_out, T_hot_in, T_hot_out):
    P = (T_cold_out-T_cold_in)/(T_hot_in-T_cold_in)
    R = (T_hot_in-T_hot_out)/(T_cold_out-T_cold_in)

    if (1-P) <= 0 or (1-P*R) <= 0 or R <= 0 or P <= 0:
        return 1e-6 #invalid range
    if abs(R - 1) < 1e-8:
        return 1e-6

    F = math.sqrt(R**2 + 1)/(R-1) * np.log10((1-P)/(1-P*R)) / (np.log10(((2/P) - 1 - R + math.sqrt(R**2 + 1))/((2/P) - 1 - R - math.sqrt(R**2 + 1))))
    return F

def F_2_2N(T_cold_in, T_cold_out, T_hot_in, T_hot_out):
    P = (T_cold_out-T_cold_in)/(T_hot_in-T_cold_in)
    R = (T_hot_in-T_hot_out)/(T_cold_out-T_cold_in)

    if (1-P) <= 0 or (1-P*R) <= 0 or R <= 0 or P <= 0:
        return 1e-6 #invalid range

    F = math.sqrt(R**2 + 1)/(2*(R-1)) * np.log10((1-P)/(1-P*R)) / (np.log10(((2/P) - 1 - R + (2/P) * math.sqrt((1-P)*(1-P*R)) + math.sqrt(R**2 + 1))/((2/P) - 1 - R + (2/P) * math.sqrt((1-P)*(1-P*R)) - math.sqrt(R**2 + 1))))
    return F

#Thermal analysis
def Thermal_analysis(Re_tube, Re_sh, N_tubes_per_pass, N_tube_passes):
    Pr = constants["Pr"]
    k_w = constants["k_w"]
    d_i = constants["d_i"]
    d_o = constants["d_o"]
    L = constants["L"]
    k_tube = constants["k_tube"]
    Nu_i = 0.023 * Re_tube**0.8 * Pr**0.3
    c = 0.2 #0.15 for square pitch, 0.2 for triangular
    Nu_o = c * Re_sh**0.6 * Pr**0.3
    h_i = Nu_i * k_w/d_i
    h_o = Nu_o * k_w/d_o
    H = (1/h_i + np.pi * d_i**2 /4 * math.log(d_o/d_i) / (2*np.pi*k_tube*L) + d_i**2/(d_o**2 * h_o))**-1
    A = N_tubes_per_pass * N_tube_passes * np.pi * d_i * L
    return H, A

def heat_balance(T_cold_out, m_dot_1, m_dot_2, T_cold_in, T_hot_in, H, A, N_shell_passes, N_tube_passes):
    Cp = constants["Cp"]
    Q1 = m_dot_1 * Cp * (T_cold_out - T_cold_in)
    T_hot_out = T_hot_in - Q1 / (m_dot_2 * Cp)

    dt_lm = delta_T_lm_counter(T_cold_in, T_cold_out, T_hot_in, T_hot_out)

    if N_shell_passes == 1 and N_tube_passes % 2 == 0:
        F = F_1_2N(T_cold_in, T_cold_out, T_hot_in, T_hot_out)
        Q_LMTD = H * A * dt_lm *F

    elif N_shell_passes == 2 and N_tube_passes > 2 and N_tube_passes % 2 == 0:
        F = F_2_2N(T_cold_in, T_cold_out, T_hot_in, T_hot_out)
        Q_LMTD = H * A * dt_lm *F

    elif N_tube_passes == N_shell_passes:
        Q_LMTD = H * A * dt_lm

    else: #Must add 3/4 shell passes
        return 1e6
    
    return Q1 - Q_LMTD  # We want this to be zero

def test_prev_year():
    T_cold_in = constants["T_cold_in"]
    T_hot_in = constants["T_hot_in"]
    Cp = constants["Cp"]

    N_tube_passes = constants["N_tube_passes"]
    N_shell_passes = constants["N_shell_passes"]
    N_rows = constants["N_rows"]
    Y = constants["Y"]

#    N_tubes_per_pass, N_B = x
#    N_B = int(round(N_B))
#    N_tubes_per_pass = int(round(N_tubes_per_pass))
    N_B = constants["N_B"]
    N_tubes_per_pass = constants["N_tubes_per_pass"]

    m_dot_1_init = constants["m_dot_1_init"]
    m_dot_2_init = constants["m_dot_2_init"]
    
    #m_dot_1_adj, m_dot_2_adj, delta_p_1, delta_p_2, pressure_cold, pressure_hot = find_flow_rates(m_dot_1_init, m_dot_2_init, int(N_tubes_per_pass), int(N_tube_passes), int(N_shell_passes), int(N_B), Y, N_rows)
    m_dot_1_adj = constants["m_dot_1"]
    m_dot_2_adj = constants["m_dot_2"]
    H, A = Thermal_analysis(Re_tube_calc(m_dot_2_adj, int(N_tubes_per_pass)), Re_sh_calc(m_dot_1_adj, int(N_B), Y, N_shell_passes), N_tubes_per_pass, N_tube_passes)

    f_a = heat_balance(T_cold_in + 1e-3, m_dot_1_adj, m_dot_2_adj, T_cold_in, T_hot_in, H, A, N_shell_passes, N_tube_passes)
    f_b = heat_balance(T_hot_in - 1e-3, m_dot_1_adj, m_dot_2_adj, T_cold_in, T_hot_in, H, A, N_shell_passes, N_tube_passes)

#    print(f_a, f_b)

    if f_a * f_b > 0:
#        print("No sign change in bracket, skipping this configuration.")
        return 1e6

    heat_calculation = root_scalar(
        heat_balance, 
        args=(m_dot_1_adj, m_dot_2_adj, T_cold_in, T_hot_in, H, A, N_shell_passes, N_tube_passes),
        bracket=[T_cold_in + 1, T_hot_in - 1],
        method='brentq',
    )

    print(f"Cold mass flow: {round(m_dot_1_adj, 4)}, Hot mass flow: {round(m_dot_2_adj, 4)}")
    #print("Pressure cold:", round(pressure_cold, 4), ", Pressure Hot:", round(pressure_hot, 4))
    print("IN FUNC --> H =", H, ", A =",A)

    if heat_calculation.converged:
        T_cold_out = heat_calculation.root
        Q = m_dot_1_adj * Cp * (T_cold_out-T_cold_in)
        T_hot_out = T_hot_in - Q / (m_dot_2_adj * Cp)
#        print(round(T_cold_out-273.15,4), round(T_hot_out-273.15, 4))
        if T_cold_out>T_hot_out:
            return 1e6 #invalid solution
        return -Q
    else:
        return 1e6

#Bounds for outlet temperatures, must lie between inlet temps
bounds = [(1, 18), (1,20)]

# if __name__ == "__main__":

#     result = differential_evolution(search, bounds, tol=1e-6)
#     print(result.x)
#     Q_dot_LMTD = -result.fun
#     print("Q_max_LMTD:", Q_dot_LMTD)
    
# ### e-NTU calcs  
# #    N_tubes_per_pass = result.x[0]
#     N_tubes_per_pass = constants["N_tubes_per_pass"]
# #    N_B = result.x[1]
#     N_B = constants["N_B"]                                  ###just to test specific years data
#     N_tube_passes = constants["N_tube_passes"]
#     N_shell_passes = constants["N_shell_passes"]
#     N_rows = constants["N_rows"]
#     Y = constants["Y"]
#     m_dot_1_init = constants["m_dot_1_init"]
#     m_dot_2_init = constants["m_dot_2_init"]
#     m_dot_1_adj, m_dot_2_adj, delta_p_1, delta_p_2, pressure_cold, pressure_hot = find_flow_rates(m_dot_1_init, m_dot_2_init, int(N_tubes_per_pass), int(N_tube_passes), int(N_shell_passes), int(N_B), Y, N_rows)
#     H, A = Thermal_analysis(Re_tube_calc(m_dot_2_adj, int(N_tubes_per_pass)), Re_sh_calc(m_dot_1_adj, int(N_B), Y, N_shell_passes), N_tubes_per_pass, N_tube_passes)
#     Q_dot_eNTU = NTU_method(m_dot_1_adj, m_dot_2_adj, H, A, N_shell_passes)
#     print("Q_dot (ENTU):", Q_dot_eNTU)
    


#square, Q_max = 13083. triangular = 14131
#1sp 1tb Q_max = 15926, 18 tubes 3 baffle
#1sp 2tb Q_max = 26584, 18 tubes 2 baffle
#2sp 1tb Q_max = 26556, 18 tubes 1 baffle
#2sp 2tb Q_max = 27005, 18 tubes 1 baffle

#correction factor F needed for multiple passes (whenever not in counterflow)
#dp_ends varies with tube passes
#d_chich_sh varies with shell passes
#add baffle losses

#reynolds dependance for baffle loss on  cold side