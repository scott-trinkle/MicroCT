import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy import interpolate
from scipy.ndimage.filters import laplace
from libtiff import TIFF


def read_data(fn):
    data = np.loadtxt(fn, dtype=np.float64)
    data = data.reshape((data.size // 2, 2))
    return data[:, 0], data[:, 1]


def log_interp(xx, yy, xx_new, kind='linear'):
    # Import model input, model output, returns log interpolation function

    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, kind=kind)

    def log_interp(zz): return np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp(xx_new)


class Material(object):

    def __init__(self, f1fn=None, f2fn=None):

        if f1fn is not None:
            self.E1, self.f1 = read_data(f1fn)

            if f2fn is None:
                self.E = self.E1

        if f2fn is not None:
            self.E2, self.f2 = read_data(f2fn)

            if self.E1.size > self.E2.size:
                self.E = self.E1
                self.f2 = log_interp(self.E2, self.f2, self.E)
            elif self.E2.size > self.E1.size:
                self.E = self.E2
                self.f1 = log_interp(self.E1, self.f1, self.E)
            elif self.E1.size == self.E2.size:
                self.E = self.E1

    def interpolate_to(self, new_E):
        self.f1_int = log_interp(self.E, self.f1, new_E)
        self.f2_int = log_interp(self.E, self.f2, new_E)
        self.E_int = new_E


class Spectrum(object):

    def __init__(self, fn):
        self.E, self.I0 = read_data(fn)
        h = 4.1357e-18  # keV s
        c = 299792458  # m / s
        self.lam = h * c / self.E  # m

    def interpolate_to(self, new_E):
        self.I0_int = log_interp(self.E, self.I0, new_E)
        self.E_int = new_E
        h = 4.1357e-18  # keV s
        c = 299792458  # m / s
        self.lam_int = h * c / new_E  # m


# get filenames of all f1, f2 and spectra data
fns = glob('../Data/Spectra/*.txt')

# Reading in energy spectra
Al = Spectrum(fns[0])
Ti = Spectrum(fns[9])

# Reading in f1 and f2 values
Os = Material(fns[4], fns[7])
U = Material(fns[5], fns[8])

# H20 only has f2 values
H2O = Material(fns[2], fns[5])

r_e = 2.818e-15  # m
R2 = 0.30  # m
pix_size = 1.24e-6  # m

Na = 6.022e23
Os.A = 190.23
U.A = 238.02891
H2O.A = 18.01528

Al_phant_tif = TIFF.open('../Data/Projection_Tifs/Al_1080.tif', 'r')
Al_phant = Al_phant_tif.read_image()

'''
NOTE

GET HISTOGRAM OF LAPLACIAN OF REAL SINOGRAM IMAGE
CONVERT INTO A "PDF" AND APPLY IT TO THE MEAN VALUE OF 
N_{A,I} FOR EACH MATERIAL TO GET A "SPECTRUM" OF 
LAPLACE-PHI(E) VALUES
'''

lap = laplace(phantom)

plt.hist(lap.flatten())
plt.show()

# mask_file = TIFF.open('../Data/Mask.tif', 'r')
# mask = mask_file.read_image() // 255

# p_water = phantom.copy() * mask  # water

# p_metal = phantom.copy()  # metal
# p_metal[mask == 1] = 0

for metal in [Os, U]:
    for spect in [Al, Ti]:

        E_min = np.array([mat.E.min() for mat in [metal, spect]]).max()
        E_max = np.array([mat.E.max() for mat in [metal, spect]]).min()
        n = np.array([mat.E.size for mat in [metal, spect]]).max()
        E = np.linspace(E_min, E_max, n)
        metal.interpolate_to(E)
        spect.interpolate_to(E)
