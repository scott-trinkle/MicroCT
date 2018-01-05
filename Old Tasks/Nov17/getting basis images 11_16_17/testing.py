'''
This script is attempting to construct U, Os and H2O basis images from the
L-edge bracketed VSO169 images
'''
import numpy as np
import matplotlib.pyplot as plt
from libtiff import TIFF
from scipy.interpolate import interp1d
from microct.phasefunctions import read_data
import numba


def du(u_func, EL, EH):
    '''
    Returns the difference in u/p between two energies

    u_func : NumPy function, takes E in keV and returns u/p in cm^2/g
    EL, EH : low and high energies for difference
    '''
    return u_func(EH) - u_func(EL)


def other_subtract(g_H, g_L, du):
    '''
    Based on decomposition of sample into a given metal and "other"
    Assumes that du for the "other" material is negligible
    '''
    return (g_H - g_L) / du


def twobytwo(g_L, g_H, u1, u2, EL, EH):
    '''
    g_L, g_H: low and high energy images
    u1, u2: u/p generating functions for material 1 and 2
    EL, EH: acquisition energies for g_L and g_H

    A[0], A[1]: Basis images for material 1 and 2
    '''
    U_mat_inv = np.linalg.inv(np.array([[u1(EL), u2(EL)],
                                        [u1(EH), u2(EH)]]))

    G = np.stack((g_L, g_H), axis=0)
    A = np.einsum('mn,nop->mop', U_mat_inv, G)
    return A[0, :, :], A[1, :, :]


@numba.jit("float64[:,:,:](float64[:,:], int16[:,:,:])")
def fourbythree_jit(U, G):
    '''
    g_L, g_H: low and high energy images2
    u1, u2: u/p generating functions for material 1 and 2
    EL, EH: acquisition energies for g_L and g_H

    A[0], A[1]: Basis images for material 1 and 2
    '''

    a = np.linalg.lstsq(U_mat, G.reshape((G.shape[0], -1)))[0]
    a = a.reshape((a.shape[0], G.shape[1], G.shape[2]))

    return a


def fourbythree_py(U, G):
    '''
    g_L, g_H: low and high energy images2
    u1, u2: u/p generating functions for material 1 and 2
    EL, EH: acquisition energies for g_L and g_H

    A[0], A[1]: Basis images for material 1 and 2
    '''

    a = np.linalg.lstsq(U_mat, G.reshape((G.shape[0], -1)))[0]
    a = a.reshape((a.shape[0], G.shape[1], G.shape[2]))

    return a


def fourbythree_test(U_mat, G):
    z, rows, cols = G.shape
    a = np.empty((z - 1, rows, cols))

    for row in range(rows):
        for col in range(cols):
            a[:, row, col] = np.linalg.lstsq(U_mat, G[:, row, col])[0]
    return a


# E in keV, u/p in cm^2 / g
NIST_path = '../../../Data/NIST/'
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


# Importing a single PROJECTION slice of the corrected data at all 4 energies
tif_path = '../l-edge raw subtraction 11_15_17/new_subs/tifs/'
G_Os_L = TIFF.open(tif_path + 'Os_low_sino_225.tif', 'r').read_image()
G_Os_H = TIFF.open(tif_path + 'Os_high_sino_225.tif', 'r').read_image()
G_U_L = TIFF.open(tif_path + 'U_low_sino_225.tif', 'r').read_image()
G_U_H = TIFF.open(tif_path + 'U_high_sino_225.tif', 'r').read_image()

n = 50
G_Os_L = G_Os_L[:n, :n]
G_Os_H = G_Os_H[:n, :n]
G_U_L = G_U_L[:n, :n]
G_U_H = G_U_H[:n, :n]

# # ASSUMING HERE THAT IMAGES ARE THE SAME SIZE
# rows, cols = G_Os_L.shape

# '''
# Four methods:
# (1) metal and "other"
# (2) 2x2 metal and water
# (3) 2x2 metal and metal
# (4) 4x3 all three
# '''

# # (1) Metal and "other"


# # Returns "basis images" based on metal and "other" decomposition
# a1U = other_subtract(G_U_H, G_U_L, du(U, EUL, EUH))
# a1Os = other_subtract(G_Os_H, G_Os_L, du(Os, EOsL, EOsH))

# print('writing v1')
# a1U.tofile('v1_U_basis')
# a1Os.tofile('v1_Os_basis')


# # (2) 2x2 Metal and "water"


# print('Calculating Os, U and water for both...')

# a2OsH2O, aOsOs = twobytwo(G_Os_L, G_Os_H, H2O, Os, EOsL, EOsH)
# a2UH2O, a2UU = twobytwo(G_U_L, G_U_H, H2O, U, EUL, EUH)

# print('writing v2')
# a2OsH2O.tofile('v2_OsIm_H2O_basis')
# aOsOs.tofile('v2_OsIm_Os_basis')
# a2UH2O.tofile('v2_UIm_H2O_basis')
# a2UU.tofile('v2_UIm_U_basis')


# # PROOF THAT THIS WORKS
# # a = G_Os_L
# # b = G_Os_H
# # u_mat = np.linalg.inv(np.array([[H2O(EOsL), Os(EOsL)],
# #                                 [H2O(EOsH), Os(EOsH)]]))

# # y, x = 1100, 3
# # G = np.array([a[y, x], b[y, x]])

# # test = np.matmul(u_mat, G)

# # print('Water: {}\n Os: {}\n'.format(aOs_wat[y, x], aOs_Os[y, x]))
# # print(test)


# # (3) 2x2 Metal and Metal

# print('\nCalculating Os, U and for both...')

# a3OsOs, a3OsU = twobytwo(G_Os_L, G_Os_H, Os, U, EOsL, EOsH)
# a3UOs, a3UU = twobytwo(G_U_L, G_U_H, Os, U, EUL, EUH)

# print('writing v3')
# a3OsOs.tofile('v3_OsIm_Os_basis')
# a3OsU.tofile('v3_OsIm_U_basis')
# a3UOs.tofile('v3_UIm_Os_basis')
# a3UU.tofile('v3_UIm_U_basis')

# # (4) 4x3 all three

# print('\n Calculating Os, U and H2O with all three images')

U_mat = np.array([[H2O(EOsL), Os(EOsL), U(EOsL)],
                  [H2O(EOsH), Os(EOsH), U(EOsH)],
                  [H2O(EUL), Os(EUL), U(EUL)],
                  [H2O(EUH), Os(EUH), U(EUH)]])


G = np.array([G_Os_L, G_Os_H, G_U_L, G_U_H])


# a = fourbythree(U_mat, G)
# print('writing v4')
# a[0].tofile('v4_H2O_basis')
# a[1].tofile('v4_Os_basis')
# a[2].tofile('v4_U_basis')
