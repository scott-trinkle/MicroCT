import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import zoom


def copync(infn, outfn):

    inf = Dataset(infn, "r")
    outf = Dataset(outfn, "w", format='NETCDF3_CLASSIC')

    #  Copy global attributes.

    outf.setncatts(inf.__dict__)

    #  Copy dimensions.

    for dimname, dim in inf.dimensions.items():
        outf.createDimension(dimname, len(dim))

    #  Copy variables and variable attributes.

    for varname, ncvar in inf.variables.items():
        var = outf.createVariable(varname, ncvar.dtype, ncvar.dimensions)
        attdict = ncvar.__dict__
        var.setncatts(attdict)
        var[:] = ncvar[:]
        outf.sync()

    inf.close()

    return outf


def downsample(infn, outfn):

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
    nz, ny, nx = 1, 5, 5
    reductions = (1 / nz, 1 / ny, 1 / nx)
    sampled = zoom(orig_data, reductions, order=1)

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


downsample('test.nc', 'new.nc')
