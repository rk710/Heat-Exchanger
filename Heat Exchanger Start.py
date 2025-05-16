import numpy as np

#Hydraulic Design
T_cold_in = 293.15
T_hot_in = 333.15
Cp = 
rho_w = 
k_w = 
mu =
Pr = 
k_tube =
L =
d_sh =
d_noz =
d_i = 0.006
d_o = 0.008
N = 
N_B =
B = N_B/L
Y =

#Hot side analysis
m_dot_1 = 
m_dot_2 = 
m_dot_tube = m_dot_2/N
v_tube = m_dot_tube*4/(rho_w*np.pi*d_i**2)
Re_tube = rho_w*v_tube*d_i/mu #check if d_i or L
v_noz_2 = m_dot_2*4/(rho_w*np.pi*d_noz**2)
f=(1.82*log(Re_tube, 10) - 1.64)**-2
delta_p_tube = 
sigma = 
K_c =
K_e =
delta_p_ends = 0.5*rho_w*v_tube**2 * (K_c + K_e)
delta_p_noz_2 = 
delta_p_2 = delta_p_tube + delta_p_ends + delta_p_noz_2

#Cold side analysis
A_sh = 
v_sh =
d_chic_sh = #Varies with number of passes
Re_sh = 
alpha = #0.2 triangular pitch, 0.34 square
delta_p_sh = 4*alpha*Re_sh**-0.15 *N*rho_w*v_sh**2
v_noz_1 = 
delta_p_noz_1 = 
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
