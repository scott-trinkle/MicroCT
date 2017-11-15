from microCTtools import get_slices
import numpy as np
from libtiff import TIFF


def subtract(fn_low, fn_high, met_name, recon, dz, n):

    print('Reading low slices')
    low_slices = get_slices(fn_low, n)
    print('Reading high slices')
    high_slices = get_slices(fn_high, n)
    print('Subtracting slices')
    subtracted = high_slices - low_slices

    slice_nums = np.arange(0, dz, dz // n)

    for i in range(n):
        print('Saving {}_low_{}'.format(met_name, slice_nums[i]))
        Al_tiff = TIFF.open(
            '../Data/new_subs/{}_low_{}_{}.tif'.format(met_name, recon, slice_nums[i]), 'w')
        Al_tiff.write_image(low_slices[i])
        Al_tiff.close()

        print('Saving {}_high_{}'.format(met_name, slice_nums[i]))
        No_tiff = TIFF.open(
            '../Data/new_subs/{}_high_{}_{}.tif'.format(met_name, recon, slice_nums[i]), 'w')
        No_tiff.write_image(high_slices[i])
        No_tiff.close()

        print('Saving {}_subtracted_{}'.format(met_name, slice_nums[i]))
        Subtraction_tiff = TIFF.open(
            '../Data/new_subs/{}_subtracted_{}_{}.tif'.format(met_name, recon, slice_nums[i]), 'w')
        Subtraction_tiff.write_image(subtracted[i])
        Subtraction_tiff.close()


# sinogram
U_low = '../Data/Brain/VS0169_1712/VS0169_1712_.volume'
U_high = '../Data/Brain/VS0169_1722/VS0169_1722_.volume'
subtract(U_low, U_high, 'U', 'sino', 900, 4)

# recon
U_low = '../Data/Brain/VS0169_1712/VS0169_1712_recon.volume'
U_high = '../Data/Brain/VS0169_1722/VS0169_1722_recon.volume'
subtract(U_low, U_high, 'U', 'recon', 1200, 4)

# sinogram
Os_low = '../Data/Brain/VS0169_Os/VS0169_Os_B_.volume'
Os_high = '../Data/Brain/VS0169_Os/VS0169_Os_A_recon.volume'
subtract(Os_low, Os_high, 'Os', 'sino', 900, 4)

# recon
Os_low = '../Data/Brain/VS0169_Os/VS0169_Os_B_recon.volume'
Os_high = '../Data/Brain/VS0169_Os/VS0169_Os_A_recon.volume'
subtract(Os_low, Os_high, 'Os', 'recon', 900, 4)
