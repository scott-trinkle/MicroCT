import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import zoom
import time

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
    
    t1 = time.time()
    print('Reading in full data...')
    orig_data = orig.variables[data_str][:]
    t2 = time.time()
    print('Done: {} seconds elapsed'.format(t2-t1))

    reductions = (1 / factors[0], 1 / factors[1], 1 / factors[2])

    print('Sampling...')
    sampled = zoom(orig_data, reductions, order=0)
    t3 = time.time()
    print('Done: {} seconds elapsed'.format(t3-t2))

    new.createDimension('dz', sampled.shape[0])
    new.createDimension('dy', sampled.shape[1])
    new.createDimension('dx', sampled.shape[2])

    new_data = new.createVariable('data', np.int16, ('dz', 'dy', 'dx'))
    attdict = orig.variables[data_str].__dict__
    new_data.setncatts(attdict)
    new_data[:] = sampled

    
    print('Writing new file...')
    new.close()
    t4 = time.time()
    print('Done: {} seconds elapsed'.format(t4-t3))
