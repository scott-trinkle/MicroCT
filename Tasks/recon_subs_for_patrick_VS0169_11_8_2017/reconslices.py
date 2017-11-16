from microCTtools import get_slices
import numpy as np
from libtiff import TIFF

fns = ['../../VS0169_Al_filterrecon.volume', '../../VS0169_No_filterrecon.volume']

n = 5
print('Reading Al slices')
Al_slices = get_slices(fns[0], n)
print('Reading No slices')
No_slices = get_slices(fns[1], n)
print('Subtracting slices')
subtracted = No_slices - Al_slices

dz = 1200
slice_nums = np.arange(0, dz, dz // n)

for i in range(n):
    print('Saving Al_{}'.format(slice_nums[i]))
    Al_tiff = TIFF.open('../../slices/Al_{}.tif'.format(slice_nums[i]), 'w')
    Al_tiff.write_image(Al_slices[i])
    Al_tiff.close()
    print('Saving No_{}'.format(slice_nums[i]))
    No_tiff = TIFF.open('../../slices/No_{}.tif'.format(slice_nums[i]), 'w')
    No_tiff.write_image(No_slices[i])
    No_tiff.close()
    print('Saving subtracted_{}'.format(slice_nums[i]))
    Subtraction_tiff = TIFF.open('../../slices/subtracted_{}.tif'.format(slice_nums[i]), 'w')
    Subtraction_tiff.write_image(subtracted[i])
    Subtraction_tiff.close()
