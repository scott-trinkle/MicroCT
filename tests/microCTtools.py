'''
This module contains functions used in processing microCT data from APS 
'''

import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import zoom
from libtiff import TIFF
import time

def downsample(infn, outfn, factors):
    '''
    This script reads in data from a netCDF file (either .nc or .volume), 
    downsamples it along all three dimensions using nearest neighbors interpolation
    according to specified factors and saves the result to a new netCDF file. 

    Parameters:
    __________
    infn : str
        Path to input, full resolution file
    outfn : str
        Path to output, downsampled file
    factors : tuple
        Downsample factors along all 3 dimensions

    Returns
    _______
    Writes sampled data to new netCDF specified by outfn

    Example
    _______
    >> downsample('data_full.nc', 'data_sampled.nc', (2, 4, 4))
    
    This will result in a new file 'data_sampled.nc' with the data from
    'data_full.nc' scaled by a factor of 2, 4 and 4 in the 1st, 2nd and
    3rd dimensions, respectively. 
    '''

    ftype_str = infn.split('.')[-1]

    if ftype_str == 'nc':
        data_str = 'array_data'
    elif ftype_str == 'volume':
        data_str = 'VOLUME'
    else:
        print('Please enter a ".nc" or a ".volume" file')
        return

    orig = Dataset(infn, 'r')
    new = Dataset(outfn, 'w', format='NETCDF3_CLASSIC')

    new.setncatts(orig.__dict__) # copies attributes
    
    t1 = time.time()
    print('Reading in full data...')
    orig_data = orig.variables[data_str][:]
    t2 = time.time()
    print('Done: {} seconds elapsed'.format(t2-t1))

    reductions = (1 / factors[0], 1 / factors[1], 1 / factors[2])

    print('Sampling...')
    sampled = zoom(orig_data, reductions, order=0) # order=0 specifies nearest neighbors
    t3 = time.time()
    print('Done: {} seconds elapsed'.format(t3-t2))

    # netCDF makes you specify variable dimensions separately...
    new.createDimension('dz', sampled.shape[0])
    new.createDimension('dy', sampled.shape[1])
    new.createDimension('dx', sampled.shape[2])

    # output netCDF files should only have one non-empty variable: "data"
    new_data = new.createVariable('data', np.int16, ('dz', 'dy', 'dx'))
    attdict = orig.variables[data_str].__dict__
    new_data.setncatts(attdict)
    new_data[:] = sampled

    
    print('Writing new file...')
    new.close()
    t4 = time.time()
    print('Done: {} seconds elapsed'.format(t4-t3))


def shift_subtract(fn, shift = (2,2)):
    '''
    This function opens a TIF file, creates a copy shifted by a specified number of pixels,
    subtracts the original image from the shifted image and returns the subtraction

    Parameters
    __________
    fn : str
        Path to input TIFF file
    shift : tuple
        Number of pixels (row, col) used to shift image in each direction. 

    Returns
    _______
    subtraction : ndarray
        Subtraction of the shifted and original image
    '''

    im_file = TIFF.open(fn, 'r')
    im = im_file.read_image()
    im_file.close()

    shape = im.shape
    new_shape = (shape[0] + shift[0], shape[1] + shift[1])

    orig_im = np.zeros(new_shape)
    shifted_im = np.zeros(new_shape)
    orig_im[:shape[0], :shape[1]] = im
    shifted_im[shift[0]:, shift[1]:] = im

    subtraction = (shifted_im - orig_im).astype(np.int16)
    return subtraction


def get_slices(fn, n):
    '''
    This function returns n evenly spaced slices from a netCDF sinogram or recon volume.

    Parameters:
    __________
    fn : str
        Path to ".volume" data
    n : int
        Number of slices to return 

    Returns:
    ________
    slices : ndarray
        Contains the n slices of the volume
    '''

    volume = nc.Dataset(fn, 'r').variables['VOLUME'][:]
    dz = volume.shape[0]
    slices = volume[0:dz:dz // 5]
    return slices
