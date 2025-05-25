import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution, root_scalar

#Compressor data (cold side)
cold_mass_flow = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0])
cold_pressure_rise = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

#Compressor data (hot side)
hot_flow_lps = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])
hot_press_bar = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")

#Hydraulic Design
constants = {

    "T_cold_in": 293.15,
    "T_hot_in": 333.15,
    "Cp": 4179,
    "rho_w": 990.1,
    "k_w": 0.632,
    "mu": 0.000651,
    "Pr": 4.31,
    "k_tube": 386,
    "d_sh": 0.064,
    "d_noz": 0.02,
    "d_i": 0.006,
    "d_o": 0.008,
    "m_dot_1_init": 0.5,
    "m_dot_2_init": 0.45,
#   "N_init": 13,
    "N_B_init": 9,
    "Y_init": 0.014,
    "N_tube_init": 1,
    "N_shell_init": 1,
    "L": 0.35,
    #### NO LONGER DEFINE N tubes, but N = N tubes/pass * N tube passes ###
    "N_tube_passes": 1,
    "N_tubes_per_pass": 13,
    "N_shell_passes": 1,
    "B_window": 1/3, #include in calcs for baffle pressure loss
    
}