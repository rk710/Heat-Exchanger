import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution

#Hydraulic Design
T_cold_in = 293.15
T_hot_in = 333.15
Cp = 4179
rho_w = 990.1
k_w = 0.632
mu = 0.000651
Pr = 4.31
k_tube = 386
L = 0.35 #Can be modified
d_sh = 0.064
d_noz = 0.02
d_i = 0.006
d_o = 0.008
N = 13 #Number of tubes, can be modified
N_B = 9 #no. of baffles
B = L/(N_B+1) #baffle spacing
Y = 0.014 #Pitch, can be modified

#Compressor data (cold side)
cold_mass_flow = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0])
cold_pressure_rise = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

#Compressor data (hot side)
hot_flow_lps = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])
hot_press_bar = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")

def heat_exchanger_pressure_drop(m_dot_1, m_dot_2):
    #Hot side analysis
    m_dot_tube = m_dot_2/N
    v_tube = m_dot_tube*4/(rho_w*np.pi*d_i**2)
    global Re_tube
    Re_tube = rho_w*v_tube*d_i/mu
    v_noz_2 = m_dot_2*4/(rho_w*np.pi*d_noz**2)
    f=(1.82*math.log(Re_tube, 10) - 1.64)**-2
    delta_p_tube = 0.5 * rho_w * v_tube**2 * f * L/d_i
    sigma = N*d_i**2/d_sh**2
    K_c = 0.4*(1-0.75*sigma) + 0.1*10000/Re_tube
    K_e = (1-2*sigma+1*sigma**2) - sigma*0.1*10000/Re_tube #Obtain K_c and K_e from figure 8 in handout, function of sigma hence are just constant. Use turbulent values.
    delta_p_ends = 0.5*rho_w*v_tube**2 * (K_c + K_e) #Must modify for multiple tube passes
    delta_p_noz_2 = rho_w * v_noz_2**2
    delta_p_2 = (delta_p_tube + delta_p_ends + delta_p_noz_2)/100000 #gives pressure in bar

    #Cold side analysis
    A_sh = d_sh/Y * (Y-d_o) * B
    v_sh = m_dot_1 / (rho_w*A_sh)
    d_chic_sh = d_sh * A_sh*4/(np.pi*d_sh**2) #Varies with number of shell passes
    global Re_sh
    Re_sh = rho_w * v_sh * d_chic_sh / mu                    
    alpha = 0.2 #0.2 triangular pitch, 0.34 square
    delta_p_sh = 4*alpha*Re_sh**-0.15 *N*rho_w*v_sh**2
    v_noz_1 = m_dot_1*4/ (rho_w*np.pi*d_noz**2)
    delta_p_noz_1 = rho_w * v_noz_1**2
    delta_p_1 = (delta_p_sh + delta_p_noz_1)/100000
    return delta_p_1, delta_p_2 #Return (delta_p_1, delta_p_2) in bar

def find_flow_rates(tol=1e-3, max_iter=1000):
    global m_dot_1
    m_dot_1 = 0.5 
    global m_dot_2
    m_dot_2 = 0.45 

    for i in range(max_iter):
        delta_p_1, delta_p_2 = heat_exchanger_pressure_drop(m_dot_1, m_dot_2)
        l_dot_1 = m_dot_1/(rho_w/1000)
        l_dot_2 = m_dot_2/(rho_w/1000)
        pressure_cold = float(cold_interp(l_dot_1))
        pressure_hot = float(hot_interp(l_dot_2))

        err_cold = abs(delta_p_1 - pressure_cold)
        err_hot = abs(delta_p_2 - pressure_hot)

        if err_cold < tol and err_hot < tol:
            print(f"Converged in {i} iterations")
            break

        if pressure_cold < delta_p_1:
            m_dot_1 -= 0.001
        else:
            m_dot_1 += 0.001

        if pressure_hot < delta_p_2:
            m_dot_2 -= 0.001
        else:
            m_dot_2 += 0.001
    print(m_dot_1, m_dot_2)
    return m_dot_1, m_dot_2, delta_p_1, delta_p_2, pressure_cold, pressure_hot

#LMTD Function
def delta_T_lm(T_cold_in, T_cold_out, T_hot_in, T_hot_out):
    dT1 = T_hot_in - T_cold_out
    dT2 = T_hot_out - T_cold_in
    #Avoid division by zero or log of zero
    if dT1 == dT2:
        return dT1
    if dT1 <= 0 or dT2 <= 0:
        return 1e-6  #Small positive value to penalize invalid regions
    return (dT1 - dT2) / np.log(dT1 / dT2)

#Objective function, negative Q to maximize it
def objective(x):
    T_cold_out, T_hot_out = x

    # Energy balance terms
    Q1 = m_dot_1 * Cp * (T_cold_out - T_cold_in)
    Q2 = m_dot_2 * Cp * (T_hot_in - T_hot_out)
    dt_lm = delta_T_lm(T_cold_in, T_cold_out, T_hot_in, T_hot_out)
    Q_LMTD = H * A * dt_lm

    #Ensure energy balance
    error1 = Q1 - Q2
    error2 = Q1 - Q_LMTD
    penalty = error1**2 + error2**2

    if penalty > 500:  #Penalty for energy imbalance
        return 1e6 + abs(penalty)**2

    #Want minimum value
    return penalty-Q_LMTD

find_flow_rates()

#Thermal analysis
Nu_i = 0.023 * Re_tube**0.8 * Pr**0.3
c = 0.2 #0.15 for square pitch, 0.2 for triangular
Nu_o = c * Re_sh**0.6 * Pr**0.3
h_i = Nu_i * k_w/d_i
h_o = Nu_o * k_w/d_o
H = (1/h_i + np.pi * d_i**2 /4 * math.log(d_o/d_i) / (2*np.pi*k_tube*L) + d_i**2/(d_o**2 * h_o))**-1
A = N * np.pi * d_i * L

#Bounds for outlet temperatures, must lie between inlet temps
#bounds = [(T_cold_in, T_hot_in), (T_cold_in, T_hot_in)]

#result = differential_evolution(objective, bounds, tol=1e-6)

#T_cold_out_optimal, T_hot_out_optimal = result.x
#Q_max = -result.fun  

#print(T_cold_out_optimal, T_hot_out_optimal, Q_max)

#square, Q_max = 13083. triangular = 14131

#correction factor F needed for multiple passes (whenever not in counterflow)
#dp_ends varies with tube passes
#d_chich_sh varies with shell passes