import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import zoom


def downsample(infn, outfn, factors):

    ftype_str = infn.split('.')[-1]

    if ftype_str == 'nc':
        data_str = 'array_data'
    elif ftype_str == 'volume':
        data_str = 'VOLUME'
    else:
        return

    orig = Dataset(infn, 'r')
    new = Dataset(outfn, 'w', format='NETCDF3_CLASSIC')

    new.setncatts(orig.__dict__)

    orig_data = orig.variables[data_str][:]
    reductions = (1 / factors[0], 1 / factors[1], 1 / factors[2])
    sampled = zoom(orig_data, reductions, order=0)

    new.createDimension('dz', sampled.shape[0])
    new.createDimension('dy', sampled.shape[1])
    new.createDimension('dx', sampled.shape[2])

    new_data = new.createVariable('data', np.int16, ('dz', 'dy', 'dx'))
    attdict = orig.variables[data_str].__dict__
    new_data.setncatts(attdict)
    new_data[:] = sampled

    new.close()

# filenames = open('filenames.txt', 'r')
# fn = ['VS0169/' + line.strip('\n') for line in filenames]


factors = (5, 3, 3)
downsample('../test.nc', '../new.nc', factors)
