import numpy as np
import math
from constants import *

#Q_max = C_min (T_hi - T_ci)
#NTU describes hea transfer across a surface
# NTU = UA/C_min, where we use 'H' instead of 'U'#
# A is heat transfer area
# we use inner diamter for heat transfer, WHY?
#Cr = C_min / C_max
#Cr = 1

#Shell and Tube > https://www.mathworks.com/help/hydro/ref/entuheattransfer.html
#[1] Holman, J. P. Heat Transfer. 9th ed. New York, NY: McGraw Hill, 2002.
#[2] Shah, R. K. and D. P. Sekulic. Fundamentals of Heat Exchanger Design. Hoboken, NJ: John Wiley & Sons, 2003.

#effectiveness for one shell pass, 2,4,6 tube passes
#e_1 = 2 / ( (1+C_rel+math.sqrt(1+(C_rel)**2)) * (1 + math.exp((-NTU)*math.sqrt(1+(C_rel)**2))) / (1 - math.exp((-NTU)*math.sqrt(1+(C_rel)**2))))
#N shell passes, 2N,4N,6N tube passes
#e = ((((1-(e_1 * C_rel))/(1-e_1))**N_shells) - 1)/((((1-(e_1 * C_rel))/(1-e_1))**N_shells) - C_rel)

def NTU_method(m_dot_1, m_dot_2, H, A, N_shells):
    Cp = constants["Cp"]
    T_hot_in = constants["T_hot_in"]
    T_cold_in = constants["T_cold_in"]
    NTU = (H * A) / Cp
    print("NTU =", round(NTU, 4))
    #print("mdot1 = ", m_dot_1, " mdot2 = ", m_dot_2)
    C_min = min(m_dot_1, m_dot_2) * Cp
    C_max = max(m_dot_1, m_dot_2) * Cp
    C_rel = C_min / C_max
    print("C_rel =", round(C_rel, 4))
    #1 shell pass, 2,4,6 tube passes
    e = 2 / ( 1+C_rel+((math.sqrt(1+(C_rel)**2)) * (1 + math.exp((-NTU)*math.sqrt(1+(C_rel)**2))) / (1 - math.exp((-NTU)*math.sqrt(1+(C_rel)**2)))))
    print("1st e:", e)
    if N_shells != 1:
        #N shell passes, 2N,4N,6N tube passes
        e = ((((1-(e * C_rel))/(1-e))**N_shells) - 1)/((((1-(e * C_rel))/(1-e))**N_shells) - C_rel)
        print("2nd e:", e)
    Q_dot = e * C_min * (T_hot_in - T_cold_in)
    return Q_dot