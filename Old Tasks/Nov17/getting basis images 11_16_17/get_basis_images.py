'''
This script is attempting to construct U, Os and H2O basis images from the
L-edge bracketed VSO169 images
'''
import numpy as np
from libtiff import TIFF
from scipy.interpolate import interp1d
from microct.phasefunctions import read_data


def du(u_func, EL, EH):
    '''
    Returns the difference in u/p between two energies

    u_func : NumPy function, takes E in keV and returns u/p in cm^2/g
    EL, EH : low and high energies for difference
    '''
    return u_func(EH) - u_func(EL)


def get_basis(g_H, g_L, du):
    '''
    Based on decomposition of sample into a given metal and "other"
    Assumes that du for the "other" material is negligible
    '''
    basis = (g_H - g_L) / du
    basis[basis < 0] = 0
    return basis


def make_G(slicenum, imtype='sino', path='../l-edge raw subtraction 11_15_17/new_subs/tifs/'):
    '''
    Imports and returns four tiff images corresponding to 4 different
    monochromatic acquisition energies.

    slicenum : str
       slice number in .tif file name
    imtype : str
       Either 'sino' or 'recon'. Default is 'sino'
    path : str
       Path to tifs

    G_* : ndarray
       Images. Returned in order of ascending energy.
    '''

    if imtype is not 'sino' and imtype is not 'recon':
        print('imtype={}. Please enter "sino" or "recon" for the imtype!'.format(imtype))

    which_im = imtype + '_' + slicenum + '.tif'

    G_Os_L = TIFF.open(path + 'Os_low_' + which_im, 'r').read_image()
    G_Os_H = TIFF.open(path + 'Os_high_' + which_im, 'r').read_image()
    G_U_L = TIFF.open(path + 'U_low_' + which_im, 'r').read_image()
    G_U_H = TIFF.open(path + 'U_high_' + which_im, 'r').read_image()

    return G_Os_L, G_Os_H, G_U_L, G_U_H


def save_basis(slicenum, imtype, savepath='results/tifs/'):
    osL, osH, uL, uH = make_G(slicenum, imtype=imtype)

    u_basis = get_basis(uH, uL, du(U, EUL, EUH))
    u_tiff = TIFF.open(savepath + 'u_{}_{}.tif'.format(imtype, slicenum), 'w')
    u_tiff.write_image(u_basis.astype(np.float32))
    u_tiff.close()

    os_basis = get_basis(osH, osL, du(Os, EOsL, EOsH))
    os_tiff = TIFF.open(
        savepath + 'os_{}_{}.tif'.format(imtype, slicenum), 'w')
    os_tiff.write_image(os_basis.astype(np.float32))
    os_tiff.close()

    return u_basis, os_basis



# E in keV, u/p in cm^2 / g
NIST_path = "../../../Data/NIST/"
H2O_E, H2O_u = read_data(NIST_path + 'u_p_H2O.txt')
U_E, U_u = read_data(NIST_path + 'u_p_U.txt')
Os_E, Os_u = read_data(NIST_path + 'u_p_Os.txt')

# Returns callable linear-interpolated functions for u/p of the materials
U = interp1d(U_E, U_u, kind='linear')
Os = interp1d(Os_E, Os_u, kind='linear')
H2O = interp1d(H2O_E, H2O_u, kind='linear')

# U_p = 19.1
# Os_p = 22.59

# # Acquisition energies
EOsL = 10.82
EOsH = 10.92
EUL = 17.12
EUH = 17.22

for recon_num, sino_num in zip([0, 300, 600, 900], [0, 225, 450, 675]):
    print('Saving Recons - {}\n Sino - {}\n'.format(recon_num, sino_num))
    save_basis(str(recon_num), 'recon')
    save_basis(str(sino_num), 'sino')
