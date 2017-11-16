import numpy as np
from scipy import interpolate
from scipy.ndimage.filters import laplace
from libtiff import TIFF


class Material(object):

    '''
    Holds properties of sample materials.
    Called with filenames to f1 and f2 data.

    Note:
    - f1 and f2 data were taken from https://physics.nist.gov/PhysRefData/FFast/html/form.html
      using the subrange 0.2-100 keV, based on the available data for the x-ray
      spectra.

    - f1 data was not available for H2O. It was created by interpolating the H
      and O data to the same E values and using:
          f1_H20 = 2 * f1_H + 1 * f1_O
      based on the Jacobsen, Kirz, Howells "X-ray Physics" chapter
    '''

    def __init__(self, f1fn=None, f2fn=None):
        '''
        Reads in f1 and f2 data. If array sizes are different, interpolates
        smaller array to the larger array, so that both can be called with the
        same energy array "E"
        '''

        self.E1, self.f1 = read_data(f1fn)
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
        '''
        Creates new interpolated f1 and f2 values for a new E vector, new_E
        '''
        self.f1_int = log_interp(self.E, self.f1, new_E)
        self.f2_int = log_interp(self.E, self.f2, new_E)


class Spectrum(object):

    '''
    Holds properties of input spectra. Called with filename to .txt spectrum
    data.

    Note:
    - Spectrum data was acquired from Figure 3 in "High-speed tomography using
      pink beam at GeoSoilEnviroCARS" by Mark Rivers, using WebPlotDigitizer
      https://automeris.io/WebPlotDigitizer/
    '''

    def __init__(self, fn):
        '''
        Reads in spectrum data and creates new attributes "lam" and "K",
        representing wavelength and wave number, respectively.
        '''

        self.E, self.I0 = read_data(fn)
        h = 4.135667662e-18  # keV s https://en.wikipedia.org/wiki/Planck_constant
        c = 299792458  # m / s https://en.wikipedia.org/wiki/Speed_of_light
        self.lam = h * c / self.E  # m
        self.K = 2 * np.pi / self.lam  # 1 / m

    def interpolate_to(self, new_E):
        self.I0_int = log_interp(self.E, self.I0, new_E)
        self.E_int = new_E
        h = 4.135667662e-18  # keV s
        c = 299792458  # m / s
        self.lam_int = h * c / new_E  # m
        self.K_int = 2 * np.pi / self.lam_int  # 1 / m


def read_data(fn):
    # Reads .txt filename, exports first (energy) and second (f1, etc)
    # columns as two numpy arrays

    data = np.loadtxt(fn, dtype=np.float64)
    data = data.reshape((data.size // 2, 2))
    return data[:, 0], data[:, 1]


def log_interp(xx, yy, xx_new, kind='linear'):
    # Import model input, model output, new input, returns log-interpolated
    # output

    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, kind=kind)

    def log_interp(zz): return np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp(xx_new)


def make_H2O_f1():
    '''
    Interpolates H f1 data (size=94) to O f1 data (size=99) and creates
    H2O f1 data.
    '''
    E_H, f1_H = read_data('../Data/Spectra/f1_H.txt')
    E_O, f1_O = read_data('../Data/Spectra/f1_O.txt')
    f1_H = log_interp(E_H, f1_H, E_O)
    f1_H2O = 2 * f1_H + 1 * f1_O
    hdr = 'Data for H20, E = 0.2 - 100 keV \n E              f1 \n keV           e atom - 1 \n INTERPOLATED FROM H AND O'
    np.savetxt('../Data/Spectra/f1_H2O_interp.txt',
               np.stack([E_O, f1_H2O], axis=1), header=hdr, comments='#')
    return


def match_arrays(H2O, spect, metal):
    # Generates a common E vector for a given spectrum, metal and water.

    # Finds the largest "minimum" energy and smallest "maximum" energy
    # for all materials, since you cannot interpolate outside of the
    # given data range.
    E_min = np.array([mat.E.min() for mat in [H2O, metal, spect]]).max()
    E_max = np.array([mat.E.max() for mat in [H2O, metal, spect]]).min()

    # Also finding the biggest "n" so that we don't lose any sharp detail
    n = np.array([mat.E.size for mat in [H2O, metal, spect]]).max()
    E = np.logspace(np.log10(E_min), np.log10(E_max), n)

    metal.interpolate_to(E)
    spect.interpolate_to(E)
    H2O.interpolate_to(E)
    return E, H2O, spect, metal


def laplace_spectrum(im, material, nbins):
    # Takes an image array, Material object and number of bins, and returns
    # lap_spec: an array of values for the laplacian of the number density
    # projection for that material and p: the "probability density" corresponding
    # to these values.

    # input is laplacian of a real microCT projection image, normalized by the
    # mean value of that image, and multiplied by the specified number density
    # projection value.
    p, bins = np.histogram(
        laplace(im / im.mean() * material.ndens), bins=nbins, density=True)

    # np.histogram returns the bin EDGES, here we are averaging to get center
    # values.
    lap_spec = np.array([np.mean([bins[i], bins[i + 1]])
                         for i in range(nbins)])
    return lap_spec, p


def calc_lap_phi(spect, H2O, metal):
    r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]

    rows, cols = H2O.phant.shape
    ebins = H2O.f1_int.size

    lap_phi = np.zeros((rows, cols, ebins))
    for i in range(rows):
        for j in range(cols):
            lap_phi[i, j, :] = r_e * spect.lam_int * \
                (H2O.f1_int * H2O.lap[i, j] + metal.f1_int * metal.lap[i, j])

    return lap_phi


def calc_T(spect, H2O, metal):
    r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]

    rows, cols = H2O.phant.shape
    ebins = H2O.f1_int.size

    T = np.zeros((rows, cols, ebins))
    for i in range(rows):
        for j in range(cols):
            T[i, j, :] = np.exp(-2 * r_e * spect.lam_int *
                                (H2O.f2_int * H2O.phant[i, j] +
                                 metal.f2_int * metal.phant[i, j]))
    return T
