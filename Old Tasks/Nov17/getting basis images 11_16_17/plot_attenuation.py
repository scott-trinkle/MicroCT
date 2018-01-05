import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from microct.phasefunctions import read_data


# E in keV, u/p in cm^2 / g
NIST_path = "/Users/scotttrinkle/GoogleDrive/Projects/MicroCT/Data/NIST/"
H2O_E, H2O_u = read_data(NIST_path + 'u_p_H2O.txt')
U_E, U_u = read_data(NIST_path + 'u_p_U.txt')
Os_E, Os_u = read_data(NIST_path + 'u_p_Os.txt')

# Returns callable linear-interpolated functions for u/p of the materials
U = interp1d(U_E, U_u, kind='linear')
Os = interp1d(Os_E, Os_u, kind='linear')
H2O = interp1d(H2O_E, H2O_u, kind='linear')


U_p = 19.1
Os_p = 22.59
# U_p = 1
# Os_p = 1

plt.close()
plt.plot(U_E, U_u * U_p, label='U')
plt.plot(Os_E, Os_u * Os_p, label='Os')
plt.legend()
plt.xlabel('Energy [keV]')
plt.ylabel(r'$\mu_i = (\mu/\rho)_i * \rho_i$')
plt.title('Linear attenuation coefficients')
plt.savefig('linatten.png', dpi=300)
