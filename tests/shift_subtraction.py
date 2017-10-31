import numpy as np
from libtiff import TIFF
import matplotlib.pyplot as plt
from glob import glob


def shift_subtract(fn, shift):
    im_file = TIFF.open(fn, 'r')
    im = im_file.read_image()
    im_file.close()

    shape = im.shape
    new_shape = (shape[0] + shift[0], shape[1] + shift[1])

    orig_im = np.zeros(new_shape)
    shifted_im = np.zeros(new_shape)
    orig_im[:shape[0], :shape[1]] = im
    shifted_im[shift[0]:, shift[1]:] = im

    return (shifted_im - orig_im).astype(np.int16)


shift = (2, 2)
fns = glob('../Subtractions/Al_*.tif')

for fn in fns:
    slice_num = fn.split('.')[-2].split('_')[-1]
    subtraction = shift_subtract(fn, shift)
    outfn = '../Subtractions/Al_{}_shifted_{}x{}.tif'.format(
        slice_num, shift[0], shift[1])
    result = TIFF.open(outfn, 'w')
    result.write_image(subtraction)
    result.close()
