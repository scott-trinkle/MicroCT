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
H2O_E, H2O_u = read_data('../../Data/NIST/u_p_H2O.txt')
U_E, U_u = read_data('../../Data/NIST/u_p_U.txt')
Os_E, Os_u = read_data('../../Data/NIST/u_p_Os.txt')

U = interp1d(U_E, U_u, kind='linear')
Os = interp1d(Os_E, Os_u, kind='linear')
H2O = interp1d(H2O_E, H2O_u, kind='linear')

EOsL = 10.82
EOsH = 10.92
EUL = 17.12
EUH = 17.22


def du(u_func, EL, EH):
    return u_func(EH) - u_func(EL)


G_1 = TIFF.open(
    '../L_edge_raw_subtraction_VS0169_11_15_2017/new_subs/tifs/Os_low_sino_225.tif', 'r').read_image()
G_2 = TIFF.open(
    '../L_edge_raw_subtraction_VS0169_11_15_2017/new_subs/tifs/Os_high_sino_225.tif', 'r').read_image()
G_3 = TIFF.open(
    '../L_edge_raw_subtraction_VS0169_11_15_2017/new_subs/tifs/U_low_sino_225.tif', 'r').read_image()
G_4 = TIFF.open(
    '../L_edge_raw_subtraction_VS0169_11_15_2017/new_subs/tifs/U_high_sino_225.tif', 'r').read_image()


'''
Four methods:
(1) metal and "other"
(2) 2x2 metal and water
(3) 2x2 metal and metal
(4) 4x3 all three
'''

# Metal and "other"


def other(g_H, g_L, du):
    return (g_H - g_L) / du


aU_o = other(G_3, G_4, du(U, EUL, EUH))
aOs_o = other(G_2, G_1, du(Os, EOsL, EOsH))


# Metal and "water"


def twobytwo(imL, imH, u1, u2, EH, EL):
    U_mat_inv = np.linalg.inv(np.array([[u1(EH), u2(EH)],
                                        [u1(EL), u2(EL)]]))
    rows, cols = imL.shape
    A = np.zeros((2, rows, cols))
    for row in np.arange(rows):
        for col in np.arange(cols):
            g = np.array([imH[row, col], imL[row, col]])
            a = np.matmul(U_mat_inv, g)
            A[:, row, col] = a
    return A[0]


print('Calculating Os and water...')
aOs_wat = twobytwo(G_1, G_2, Os, H2O, EOsL, EOsH)
print('Calculating U and water...')
aU_wat = twobytwo(G_3, G_4, U, H2O, EUL, EUH)
