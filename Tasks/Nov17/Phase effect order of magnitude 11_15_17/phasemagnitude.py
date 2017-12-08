import numpy as np
import matplotlib.pyplot as plt
from microct.phasefunctions import *
from libtiff import TIFF
from scipy.ndimage.filters import laplace
from scipy.ndimage import zoom

specpath = '../../../Data/Spectra/'
nistpath = '../../../Data/NIST/'

# Reading in energy spectra
Al = Spectrum(specpath + 'Al_1mm_3mrad.txt')
Ti = Spectrum(specpath + 'Ti_1mm_2mrad.txt')

# Reading in f1 and f2 values
Os = Material(nistpath + 'f1_Os.txt', nistpath + 'f2_Os.txt')
U = Material(nistpath + 'f1_U.txt', nistpath + 'f2_U.txt')
H2O = Material(nistpath + 'f1_H2O_interp.txt', nistpath + 'f2_H20.txt')

# Constants
r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]
R2 = 0.30  # m [.setup files]

Na = 6.02214e23  # atoms/mole, https: // en.wikipedia.org / wiki / Avogadro_constant

Os.A = 190.23  # g / mol    https://en.wikipedia.org/wiki/Osmium
U.A = 238.02891  # g / mol  https://en.wikipedia.org/wiki/Uranium
H2O.A = 18.01488  # g / mol https://en.wikipedia.org/wiki/Molar_mass

pix_size = 1.24e-6  # m [.setup files]
width = 1920 * pix_size  # m

# Calculating number densities with formular from [Jacobsen, Kirz, Howells]
# Densities from the above cited wiki pages.
# Note these are actually the number density LINE INTEGRALS in
# particles / cm^2.
H2O.ndens = 1 * Na / H2O.A * width * 100e3  # particles / m^2
Os.ndens = 22.59 * Na / Os.A * width * 100e3  # particles / m^2
U.ndens = 19.1 * Na / U.A * width * 100e3  # particles / m^2

# Reading in image files to calculate "Laplacian Spectrums"
Al_phant = TIFF.open('../Data/Projection_Tifs/Al_1080.tif', 'r').read_image()
No_phant = TIFF.open('../Data/Projection_Tifs/No_1080.tif', 'r').read_image()

# Al_phant = zoom(Al_phant, (1 / 10, 1 / 10), order=0)
# No_phant = zoom(No_phant, (1 / 10, 1 / 10), order=0)

H2O.phant = Al_phant / Al_phant.mean() * H2O.ndens
Os.phant = No_phant / No_phant.mean() * Os.ndens
U.phant = No_phant / No_phant.mean() * U.ndens

print('\nCalculating laplacians:')
print('H2O...')
H2O.lap = laplace(H2O.phant)
print('Os...')
Os.lap = laplace(Os.phant)
print('U...\n')
U.lap = laplace(U.phant)


def forwardmodel(metal=Os, spect=Al, H2O=H2O):
    '''
    Matches all array sizes, calculates laplacian of phi and transmission factor,
    Integrates over energy to calculate the full measured intensity, Im_full
    and the transmission-only intensity Im_trans
    '''

    if metal is Os:
        nmet = 'Os'
    elif metal is U:
        nmet = 'U'
    if spect is Al:
        nspect = 'Al'
    elif spect is Ti:
        nspect = 'Ti'

    E, H2O, spect, metal = match_arrays(H2O, spect, metal)

    print('Calculating lap_phi_E...')
    lap_phi_E = calc_lap_phi(spect, H2O, metal)

    print('Calculating T_E...')
    T_E = calc_T(spect, H2O, metal)
    plt.hist(T_E.flatten(), bins=100)
    plt.show()

    return T_E

    print('Integrating Im...')
    Im_full = np.trapz(np.array([E[i] * spect.I0_int[i] *
                                 T_E[:, :, i] * (1 + R2 / spect.K_int[i] *
                                                 lap_phi_E[:, :, i]) for i in range(E.size)]), E, axis=0)

    print('Integrating Im trans...')
    Im_trans = np.trapz(np.array([E[i] * spect.I0_int[i] * T_E[:, :, i]
                                  for i in range(E.size)]), E, axis=0)

    # print('Calculating Im_phase...')
    # Im_phase = np.array([E[i] * spect.I0_int[i] *
    #                      (1 + R2 / spect.K_int[i] *
    # lap_phi_E[:, :, i]) for i in range(E.size)]).sum(axis=0)

    print('Calculating ph_factor...')
    ph_factor = np.array([(R2 / spect.K_int[i] * lap_phi_E[:, :, i])
                          for i in range(E.size)])

    return Im_full, Im_trans, ph_factor


def saveresults(result, fn):
    '''
    Saves a histogram of the results
    '''
    plt.close()
    plt.hist(result.flatten(), bins=500)
    if fn is not None:
        plt.title(fn)
        plt.xlabel('Percent Difference')
        plt.savefig('data/' +
                    fn + '.png', dpi=400)
        plt.close()
    else:
        plt.show()


def print_stats(a):
    a = a.flatten()
    print('Mean: {} \n abs Mean: {} \n Min: {} \n Max: {} \n'.format(a.mean(),
                                                                     abs(a).mean(), a.min(), a.max()))


'''
Subtraction histogram plots
'''
for metal in [Os, U]:
    '''
    For each metal, saves histogram of the percent difference between the full-
    and transmission-only measured intensities for each filter, as well as the
    percent difference in full- and transmission-only subtractions between the
    two filters. Also saves .png file of full subtraction
    '''

    plt.close()

    if metal is Os:
        nmet = 'Os'
    elif metal is U:
        nmet = 'U'

    print('{}:'.format(nmet))

    imAl, imAlT, imAlPh = forwardmodel(metal, Al)
    imTi, imTiT, imTiPh = forwardmodel(metal, Ti)

    print('Saving phase histogram...\n')
    saveresults(imAlPh, '{}/{}_phase'.format(nmet, 'Al'))
    saveresults(imTiPh, '{}/{}_phase'.format(nmet, 'Ti'))

    print('imAlPh:')
    print_stats(imAlPh)
    print('imTiPh:')
    print_stats(imTiPh)

    # print('Saving histograms...')
    # saveresults((imAl - imAlT) / imAl * 100,
    #             '{}/{}_full_vs_T'.format(nmet, 'Al'))
    # saveresults((imTi - imTiT) / imTi * 100,
    #             '{}/{}_full_vs_T'.format(nmet, 'Ti'))

    print('Energy subtraction histograms...')

    sub = imAl - imTi
    subT = imAlT - imTiT

    E_effect = (sub - subT) / sub * 100

    saveresults(E_effect, '{}/E_sub_full_vs_T'.format(nmet))
    print_stats(E_effect)

    # print('Saving image...\n')
    # plt.imshow(sub, cmap=plt.cm.Greys_r)
    # plt.axis('off')
    # plt.tight_layout()
    # plt.savefig('../Data/Spectra/Results/{}/sub_im.png'.format(nmet),
    #             dpi=1000, bbox_inches='tight')
    plt.close()
