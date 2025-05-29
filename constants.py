import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution, root_scalar

#Hydraulic Design
constants = {

    "T_cold_in": 23.7 + 273.15,
    "T_hot_in": 55.5 + 273.15,
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
    "N_B_init": 9,
    "Y": 0.012,
    "N_tube_init": 1,
    "N_shell_init": 1,
    "L": 0.25,
    #### NO LONGER DEFINE N tubes, but N = N tubes/pass * N tube passes ###
    "N_B": 8,
    "N_tube_passes": 2,
    "N_tubes_per_pass": 5,
    "N_shell_passes": 1,
    "N_rows": 4,
    "B_window": 0.25, #include in calcs for baffle pressure loss
    "Year": 2023,
    "m_dot_1": 0.6,
    "m_dot_2": 0.382,
    
}

Year = constants["Year"]
if Year == 2025:
    #Compressor data (cold side)
    cold_mass_flow = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0])
    cold_pressure_rise = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
    cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

    #Compressor data (hot side)
    hot_flow_lps = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])
    hot_press_bar = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
    hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")

elif Year == 2024:
    #Compressor data (cold side)
    cold_mass_flow = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
    cold_pressure_rise = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776])
    cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

    #Compressor data (hot side)
    hot_flow_lps = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694])
    hot_press_bar = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615])
    hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")

elif Year == 2023:
    #Compressor data (cold side)
    cold_mass_flow = np.array([0.7083, 0.6417, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
    cold_pressure_rise = np.array([0.1310, 0.2017, 0.2750, 0.3417, 0.4038, 0.4503, 0.4856, 0.5352, 0.5717, 0.5876])
    cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

    #Compressor data (hot side)
    hot_flow_lps = np.array([0.4722, 0.4340, 0.3924, 0.3507, 0.3021, 0.2523, 0.1979, 0.1493, 0.1111, 0.0694])
    hot_press_bar = np.array([0.0538, 0.1192, 0.1727, 0.2270, 0.2814, 0.3366, 0.3907, 4456, 0.4791, 0.5515])
    hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")

elif Year == 2022:
    #Compressor data (cold side)
    cold_mass_flow = np.array([0.5833, 0.5083, 0.4750, 0.4250, 0.3792, 0.3417, 0.2958, 0.2583, 0.2125, 0.1708])
    cold_pressure_rise = np.array([0.1113, 0.2157, 0.2538, 0.3168, 0.3613, 0.4031, 0.4511, 0.4846, 0.5181, 0.5573])
    cold_interp = interp1d(cold_mass_flow, cold_pressure_rise, fill_value="extrapolate")

    #Compressor data (hot side)
    hot_flow_lps = np.array([0.4583, 0.4236, 0.4010, 0.3611, 0.3125, 0.2639, 0.2222, 0.1597, 0.1181, 0.0694])
    hot_press_bar = np.array([0.1333, 0.1756, 0.2024, 0.2577, 0.3171, 0.3633, 0.4233, 0.4784, 0.5330, 0.5715])
    hot_interp = interp1d(hot_flow_lps, hot_press_bar, fill_value="extrapolate")