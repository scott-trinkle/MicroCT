'''
This script is attempting to construct U, Os and H2O basis images from the
L-edge bracketed VSO169 images
'''

import numpy as np
import matplotlib.pyplot as plt
from libtiff import TIFF
from scipy.interpolate import interp1d
from microct.phasefunctions import read_data

# E in keV, u/p in cm^2 / g
NIST_path = '/Users/scotttrinkle/GoogleDrive/Projects/MicroCT/Data/NIST/'
H2O_E, H2O_u = read_data(NIST_path + 'u_p_H2O.txt')
U_E, U_u = read_data(NIST_path + 'u_p_U.txt')
Os_E, Os_u = read_data(NIST_path + 'u_p_Os.txt')

# Returns callable linear-interpolated functions for u/p of the materials
U = interp1d(U_E, U_u, kind='linear')
Os = interp1d(Os_E, Os_u, kind='linear')
H2O = interp1d(H2O_E, H2O_u, kind='linear')

# Acquisition energies
EOsL = 10.82
EOsH = 10.92
EUL = 17.12
EUH = 17.22


def du(u_func, EL, EH):
    '''
    Returns the difference in u/p between two energies
    '''

    return u_func(EH) - u_func(EL)


# Importing a single PROJECTION slice of the corrected data at all 4 energies
tif_path = '/Users/scotttrinkle/GoogleDrive/Projects/MicroCT/Tasks/11/l-edge raw subtraction 11_15_17/new_subs/tifs/'
G_Os_L = TIFF.open(tif_path + 'Os_low_sino_225.tif', 'r').read_image()
G_Os_H = TIFF.open(tif_path + 'Os_high_sino_225.tif', 'r').read_image()
G_U_L = TIFF.open(tif_path + 'U_low_sino_225.tif', 'r').read_image()
G_U_H = TIFF.open(tif_path + 'U_high_sino_225.tif', 'r').read_image()

# ASSUMING HERE THAT IMGAES ARE THE SAME SIZE
rows, cols = G_Os_L.shape

'''
Four methods:
(1) metal and "other"
(2) 2x2 metal and water
(3) 2x2 metal and metal
(4) 4x3 all three
'''

# (1) Metal and "other"


def other_subtract(g_H, g_L, du):
    '''
    Based on decomposition of sample into a given metal and "other"
    Assumes that du for the "other" material is negligible
    '''
    return (g_H - g_L) / du


# Returns "basis images" based on metal and "other" decomposition
aU_o = other_subtract(G_U_H, G_U_L, du(U, EUL, EUH))
aOs_o = other_subtract(G_Os_H, G_Os_L, du(Os, EOsL, EOsH))


# (2) 2x2 Metal and "water"


def twobytwo(g_L, g_H, u1, u2, EL, EH):
    '''
    g_L, g_H: low and high energy images
    u1, u2: callable u/p functions for materials 1 and 2
    EL, EH: low and high energies


    '''

    # Inverse of u/p matrix
    U_mat_inv = np.linalg.inv(np.array([[u1(EL), u2(EL)],
                                        [u1(EH), u2(EH)]]))

    g = np.array([g_L, g_H])
    a = np.matmul(U_mat_inv, g)  # a[0] = u1, a[1] = u2
    return a[0], a[1]


print('Calculating Os, U and water for both...')

aOs_wat, aOs_Os = np.zeros((rows, cols)), np.zeros((rows, cols))
aU_wat, aU_U = np.zeros((rows, cols)), np.zeros((rows, cols))

time_count = 0

for row in range(rows):
    if row / rows * 100 % 10 == 0:
        time_count += 1
        print('{}% done...'.format(round(time_count / 10 * 100, 2)))
    for col in range(cols):
        aOs_wat[row, col], aOs_Os[row, col] = twobytwo(
            G_Os_L[row, col], G_Os_H[row, col], H2O, Os, EOsL, EOsH)

        aU_wat[row, col], aU_U[row, col] = twobytwo(
            G_U_L[row, col], G_U_H[row, col], H2O, U, EUL, EUH)


# (3) 2x2 Metal and Metal

print('\nCalculating Os, U and for both...')

aOs_Os, aOs_U = np.zeros((rows, cols)), np.zeros((rows, cols))
aU_Os, aU_U = np.zeros((rows, cols)), np.zeros((rows, cols))

time_count = 0

for row in range(rows):
    if row / rows * 100 % 10 == 0:
        time_count += 1
        print('{}% done...'.format(round(time_count / 10 * 100, 2)))
    for col in range(cols):
        aOs_Os[row, col], aOs_U[row, col] = twobytwo(
            G_Os_L[row, col], G_Os_H[row, col], Os, U, EOsL, EOsH)

        aU_Os[row, col], aU_U[row, col] = twobytwo(
            G_U_L[row, col], G_U_H[row, col], Os, U, EUL, EUH)
