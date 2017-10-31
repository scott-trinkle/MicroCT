import netCDF4 as nc
import numpy as np
from libtiff import TIFF


def get_sino_slices(fn):
    volume = nc.Dataset(fn, 'r').variables['data'][:]
    numslices = volume.shape[0]
    slices = volume[0:numslices:numslices // 5]
    return slices


# fns = ['../VS0169/VS0169_Al_filter.volume',
#        '../VS0169/VS0169_No_filter.volume']

fns = ['../Sampled_Data/sampled_133_VS0169_Al_filter.volume',
       '../Sampled_Data/sampled_133_VS0169_No_filter.volume']

Al_slices = get_sino_slices(fns[0])
No_slices = get_sino_slices(fns[1])

subtracted = No_slices - Al_slices

slices = np.arange(0, 1800, 1800 // 5)

for i in range(subtracted.shape[0]):
    tiff = TIFF.open('subtracted_{}.tif'.format(slices[i]), 'w')
    tiff.write_image(subtracted[i])
    tiff.close()
