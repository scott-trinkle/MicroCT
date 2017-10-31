'''
This script reads in projection TIF images from the aluminum filtered dataset and
generates the shifted subtraction images
'''

import numpy as np
from libtiff import TIFF
from glob import glob
from microCTtools import shift_subtract

shift = (2, 2)
fns = glob('../Data/Projection_Tifs/Al_*.tif')

for fn in fns:
    slice_num = fn.split('.')[-2].split('_')[-1]
    subtraction = shift_subtract(fn, shift)
    outfn = '../Data/Projection_Tifs/Al_{}_shifted_{}x{}.tif'.format(
        slice_num, shift[0], shift[1])
    result = TIFF.open(outfn, 'w')
    result.write_image(subtraction)
    result.close()
