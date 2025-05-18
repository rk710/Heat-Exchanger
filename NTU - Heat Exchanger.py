import numpy as np
import math

#Hydraulic Design
T_cold_in = 293.15
T_hot_in = 333.15
#Cp = 4179
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

#Hot side analysis
m_dot_1 = 0.5 #Initial guess
m_dot_2 = 0.45 #Initial guess
m_dot_tube = m_dot_2/N
v_tube = m_dot_tube*4/(rho_w*np.pi*d_i**2)
Re_tube = rho_w*v_tube*d_i/mu
v_noz_2 = m_dot_2*4/(rho_w*np.pi*d_noz**2)
f=(1.82*math.log(Re_tube, 10) - 1.64)**-2
delta_p_tube = 0.5 * rho_w * v_tube**2 * f * L/d_i
sigma = N*d_i**2/d_sh**2
K_c = 0.45
K_e = 0.8 #Obtain K_c and K_e from figure 8 in handout, function of sigma hence are just constant
delta_p_ends = 0.5*rho_w*v_tube**2 * (K_c + K_e)
delta_p_noz_2 = rho_w * v_noz_2**2
delta_p_2 = delta_p_tube + delta_p_ends + delta_p_noz_2

#Cold side analysis
A_sh = d_sh/Y * (Y-d_o) * B
v_sh = m_dot_1 / (rho_w*A_sh)
d_chic_sh = d_sh * A_sh*4/(np.pi*d_sh**2) #Varies with number of passes
Re_sh = rho_w * v_sh * d_chic_sh / mu
alpha = 0.34 #0.2 triangular pitch, 0.34 square
delta_p_sh = 4*alpha*Re_sh**-0.15 *N*rho_w*v_sh**2
v_noz_1 = m_dot_1*4/ (rho_w*np.pi*d_noz**2)
delta_p_noz_1 = rho_w * v_noz_1**2
delta_p_1 = delta_p_sh + delta_p_noz_1
#Need iterative process to find optimal mass flow rate

#Thermal analysis
Nu_i = 0.023 * Re_tube**0.8 * Pr**0.3
c = 0.15 #0.15 for square pitch, 0.2 for triangular
Nu_o = c * Re_sh**0.6 * Pr**0.3
#Reynolds no. adjustment?
h_i = Nu_i * k_w/d_i
h_o = Nu_o * k_w/d_o
H = (1/h_i + np.pi * d_i**2 /4 * math.log(d_o/d_i) / (2*np.pi*k_tube*L) + d_i**2/(d_o**2 * h_o))**-1
#print(B, v_tube, Re_tube, v_noz_2, delta_p_tube, delta_p_ends, delta_p_noz_2, delta_p_2, A_sh, v_sh, d_chic_sh, Re_sh, delta_p_sh, v_noz_1, delta_p_noz_1, delta_p_1, Nu_i, Nu_o, h_i, h_o, H)
#Need iteration to get T_1_out and T_2_out, correction factor F needed for multiple passes

#Q_max = C_min (T_hi - T_ci)
#NTU describes hea transfer across a surface
# NTU = UA/C_min, where we use 'H' instead of 'U'#
# A is heat transfer area
# we use inner diamter for heat transfer, WHY?
#Cr = C_min / C_max
#Cr = 1
Cp = 4179
Q_dot_max = Cp * (T_hot_in - T_cold_in)
A_tube_inner = N * np.pi * d_i * L
NTU = (H * A_tube_inner) / Cp
C_max = 4185 #at 60 degrees celsius
C_min = Cp #at midpoint of 40 degrees celsius
C_rel = C_max / C_min
N_shells = 1 #N shell passes

#Shell and Tube > https://www.mathworks.com/help/hydro/ref/entuheattransfer.html

#effectiveness for one shell pass, 2,4,6 tube passes
e_1 = 2 / ( (1+C_rel+math.sqrt(1+(C_rel)**2)) * (1 + math.exp((-NTU)*math.sqrt(1+(C_rel)**2))) / (1 - math.exp((-NTU)*math.sqrt(1+(C_rel)**2))))
#  e = ((((1-(e_1 * C_rel))/(1-e_1))**N_shells) - 1)/((((1-(e_1 * C_rel))/(1-e_1))**N_shells) - C_rel)

Q_dot = e * Q_dot_max
T_hot_out = T_hot_in - (Q_dot / Cp)
T_cold_out = T_cold_in + (Q_dot / Cp)

print(T_hot_out - 273.15)
print(T_cold_out - 273.15)