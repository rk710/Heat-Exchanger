import numpy as np
from scipy.interpolate import interp1d

# Compressor data (cold side)
cold_mass_flow = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0.0])
cold_pressure_rise = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077])
cold_interp = interp1d(cold_mass_flow, cold_pressure_rise)

# Compressor data (hot side)
hot_flow_lps = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010])
hot_press_bar = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239])
hot_interp = interp1d(hot_flow_lps, hot_press_bar)

# Function you must define
def heat_exchanger_pressure_drop(m_dot_1, m_dot_2):
    # Implement your model here
    # Return (dp_cold, dp_hot) in bar
    pass

# Iterative solver
def find_flow_rates(tol=1e-3, max_iter=100):
    m_dot_1 = 0.45 
    m_dot_2 = 0.5 
    l_dot_1 = m_dot_1/(rho_w/1000)
    l_dot_2 = m_dot_2/(rho_w/1000)

    for i in range(max_iter):
        delta_p_1, delta_p_2 = heat_exchanger_pressure_drop(m_dot_1, m_dot_2)
        pressure_cold = float(cold_interp(l_dot_1))
        pressure_hot = float(hot_interp(l_dot_2))

        err_cold = abs(delta_p_1 - pr_cold)
        err_hot = abs(delta_p_2 - pr_hot)

        if err_cold < tol and err_hot < tol:
            print(f"Converged in {i} iterations")
            break

        # Simple correction strategy: move flow in the direction that balances pressure drop
        if pressure_cold < delta_p_1:
            m_dot_1 -= 0.01
        else:
            m_dot_1 += 0.01

        if pressure_hot < delta_p_2:
            m_dot_2 -= 0.01
        else:
            m_dot_2 += 0.01

    return m_dot_1, m_dot_2, delta_p_1, delta_p_2, pressure_cold, pressure_hot
