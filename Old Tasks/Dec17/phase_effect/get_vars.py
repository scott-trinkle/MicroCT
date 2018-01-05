import numpy as np
from microct.phasefunctions import *
from libtiff import TIFF
from scipy.ndimage.filters import laplace
from scipy.ndimage import zoom

# Filepaths to data
specpath = '../../../Data/Spectra/'
nistpath = '../../../Data/NIST/'
impath = '../../Nov17/getting basis images 11_16_17/results/tifs/'

# Reading in energy spectra
Al = Spectrum(specpath + 'Al_1mm_3mrad.txt')
Ti = Spectrum(specpath + 'Ti_1mm_2mrad.txt')

# Reading in f1 and f2 values
Os = Material(nistpath + 'f1_Os.txt', nistpath + 'f2_Os.txt')
U = Material(nistpath + 'f1_U.txt', nistpath + 'f2_U.txt')

# Reading in basis images, pixels are density projections [g/cm2]
Os_phant = TIFF.open('crop_os.tif', 'r').read_image()
U_phant = TIFF.open('crop_u.tif', 'r').read_image()

# Downsample images for testing
# scale = 1
# Os_phant = zoom(Os_phant, (1 / scale, 1 / scale), order=0)
# U_phant = zoom(U_phant, (1 / scale, 1 / scale), order=0)

# Constants
r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]
R2 = 0.30  # m [.setup files]
Na = 6.02214e23  # atoms/mole, https: // en.wikipedia.org / wiki / Avogadro_constant
Os.A = 190.23  # g / mol    https://en.wikipedia.org/wiki/Osmium
U.A = 238.02891  # g / mol  https://en.wikipedia.org/wiki/Uranium

# Calculating number density projections with formula from [Jacobsen,
# Kirz, Howells]
Os.ndens = Os_phant * Na / Os.A * 100e2  # particles / m^2
U.ndens = U_phant * Na / U.A * 100e2  # particles / m^2

# Calculate laplacian of images
print('\nCalculating laplacians:')
print('Os...')
Os.lap = laplace(Os.ndens)
print('U...\n')
U.lap = laplace(U.ndens)
