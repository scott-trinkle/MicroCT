import numpy as np
import matplotlib.pyplot as plt
from phasefunctions import Material, Spectrum, laplace_spectrum, match_arrays, calc_lap_phi, calc_T
from libtiff import TIFF
from scipy.ndimage.filters import laplace
from scipy.ndimage import zoom


# Reading in energy spectra
Al = Spectrum('../Data/Spectra/Al_1mm_3mrad.txt')
Ti = Spectrum('../Data/Spectra/Ti_1mm_2mrad.txt')

# Reading in f1 and f2 values
Os = Material('../Data/Spectra/f1_Os.txt', '../Data/Spectra/f2_Os.txt')
U = Material('../Data/Spectra/f1_U.txt', '../Data/Spectra/f2_U.txt')
H2O = Material('../Data/Spectra/f1_H2O_interp.txt',
               '../Data/Spectra/f2_H20.txt')

# Constants
r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]
R2 = 0.30  # m [.setup files]

Na = 6.02214e23  # atoms/mole, https://en.wikipedia.org/wiki/Avogadro_constant

Os.A = 190.23  # g / mol    https://en.wikipedia.org/wiki/Osmium
U.A = 238.02891  # g / mol  https://en.wikipedia.org/wiki/Uranium
H2O.A = 18.01488  # g / mol https://en.wikipedia.org/wiki/Molar_mass

pix_size = 1.24e-6  # m [.setup files]
width = 1920 * pix_size * 100  # cm

# Calculating number densities with formular from [Jacobsen, Kirz, Howells]
# Densities from the above cited wiki pages.
# Note these are actually the number density LINE INTEGRALS in
# particles / cm^2.
H2O.ndens = 1 * Na / H2O.A * width
Os.ndens = 22.59 * Na / Os.A * width
U.ndens = 19.1 * Na / U.A * width

# Reading in image files to calculate "Laplacian Spectrums"
Al_phant = TIFF.open('../Data/Projection_Tifs/Al_1080.tif', 'r').read_image()
No_phant = TIFF.open('../Data/Projection_Tifs/No_1080.tif', 'r').read_image()

Al_phant = zoom(Al_phant, (1 / 10, 1 / 10), order=0)
No_phant = zoom(No_phant, (1 / 10, 1 / 10), order=0)

H2O.phant = Al_phant / Al_phant.mean() * H2O.ndens
Os.phant = No_phant / No_phant.mean() * Os.ndens
U.phant = No_phant / No_phant.mean() * U.ndens

print('')
print('Calculating laplacians:')
print('H2O...')
H2O.lap = laplace(H2O.phant)
print('Os...')
Os.lap = laplace(Os.phant)
print('U...\n')
U.lap = laplace(U.phant)


def forwardmodel(metal=Os, spect=Al, H2O=H2O):
    E, H2O, spect, metal = match_arrays(H2O, spect, metal)

    print('Calculating lap_phi_E...')
    lap_phi_E = calc_lap_phi(spect, H2O, metal)
    print('Calculating T_E...')
    T_E = calc_T(spect, H2O, metal)

    print('Calculating Im...')
    Im_full = np.array([E[i] * spect.I0_int[i] *
                        T_E[:, :, i] * (1 + R2 / spect.K_int[i] *
                                        lap_phi_E[:, :, i]) for i in range(E.size)]).sum(axis=0)
    Im_trans = np.array([E[i] * spect.I0_int[i] * T_E[:, :, i]
                         for i in range(E.size)]).sum(axis=0)

    return Im_full, Im_trans


def saveresults(result, fn):
    plt.close()
    plt.hist(result.flatten(), bins=500)
    plt.title(fn)
    plt.xlabel('Percent Difference')
    plt.savefig('../Data/Spectra/Results/' + fn + '.png', dpi=400)
    plt.close()


for metal in [Os, U]:

    if metal is Os:
        nmet = 'Os'
    elif metal is U:
        nmet = 'U'

    print('{}:'.format(nmet))

    imAl, imAlT = forwardmodel(metal, Al)
    imTi, imTiT = forwardmodel(metal, Ti)

    print('Saving histograms...')
    saveresults((imAl - imAlT) / imAl * 100,
                '{}/{}_full_vs_T'.format(nmet, 'Al'))
    saveresults((imTi - imTiT) / imTi * 100,
                '{}/{}_full_vs_T'.format(nmet, 'Ti'))

    sub = imAl - imTi
    subT = imAlT - imTiT

    saveresults((sub - subT) / sub * 100, '{}/E_sub_full_vs_T'.format(nmet))

    print('Saving image...\n')
    plt.imshow(sub, cmap=plt.cm.Greys_r)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('../Data/Spectra/Results/{}/sub_im.png'.format(nmet),
                dpi=1000, bbox_inches='tight')
    plt.close()
