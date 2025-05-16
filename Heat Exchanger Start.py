import numpy as np
import math

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

#Hot side analysis
m_dot_1 =  #Initial guess
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
Nu_i = 
Nu_o = 
#Reynolds no. adjustment?
h_i =
h_o =
H = ()**-1
#Need iteration to get T_1_out and T_2_out, correction factor F needed for multiple passes
