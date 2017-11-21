'''
This script generates raw subtraction TIFs of the Al_filter and No_filter sinograms.
'''

import numpy as np
from libtiff import TIFF
from microCTtools import get_slices

fns = ['../Data/VS0169_Sampled/sampled_133_VS0169_Al_filter.volume',
       '../Data/VS0169_Sampled/sampled_133_VS0169_No_filter.volume']

n = 5
Al_slices = get_slices(fns[0], n)
No_slices = get_slices(fns[1], n)
subtracted = No_slices - Al_slices

dz = Al_slices.shape[0] # 1800 for .volume sinograms
slices = np.arange(0, dz, dz // n)

for i in range(subtracted.shape[0]):
    tiff = TIFF.open('../Data/Projection_Tifs/subtracted_{}.tif'.format(slices[i]), 'w')
    tiff.write_image(subtracted[i])
    tiff.close()
