import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../microct/')
import phasefunctions as pf
from libtiff import TIFF
from get_vars import *


def calc_intensity(spect, Os=Os, U=U):
    '''
    Matches all array sizes, calculates laplacian of phi and transmission factor,
    Integrates over energy to calculate the full measured intensity, Im_full
    and the transmission-only intensity Im_trans
    '''

    E, spect, Os, U = pf.match_arrays(spect, Os, U)

    print('Calculating lap_phi_E...')
    lap_phi_E = pf.calc_lap_phi(spect, Os, U)

    print('Calculating T_E...')
    T_E = pf.calc_T(spect, Os, U)

    print('Integrating Im...')
    Im_full = np.trapz(np.array([E[i] * spect.I0_int[i] *
                                 T_E[i, :, :] * (1 + R2 / spect.K_int[i] *
                                                 lap_phi_E[i, :, :]) for i in range(E.size)]), E, axis=0)

    print('Integrating Im trans...')
    Im_trans = np.trapz(np.array([E[i] * spect.I0_int[i] * T_E[i, :, :]
                                  for i in range(E.size)]), E, axis=0)

    # print('Calculating Im_phase...')
    # Im_phase = np.array([E[i] * spect.I0_int[i] *
    #                      (1 + R2 / spect.K_int[i] *
    # lap_phi_E[:, :, i]) for i in range(E.size)]).sum(axis=0)

    # print('Calculating ph_factor...')
    # ph_factor = np.array([(R2 / spect.K_int[i] * lap_phi_E[:, :, i])
    #                       for i in range(E.size)])

    return Im_full, Im_trans


print('matching arrays...')
E, Al, Os, U = pf.match_arrays(Al, Os, U)

print('calculating T...')
T = pf.calc_T(Al, Os, U)

# es, row, col = T.shape

# plt.plot(E, T[:, row // 2, col // 2])
# plt.show()


# Al_full, Al_trans = forwardmodel(Al)
# Ti_full, Ti_trans = forwardmodel(Ti)


# def saveresults(result, fn):
#     '''
#     Saves a histogram of the results
#     '''
#     plt.close()
#     plt.hist(result.flatten(), bins=500)
#     if fn is not None:
#         plt.title(fn)
#         plt.xlabel('Percent Difference')
#         plt.savefig('data/' +
#                     fn + '.png', dpi=400)
#         plt.close()
#     else:
#         plt.show()


# def print_stats(a):
#     a = a.flatten()
#     print('Mean: {} \n abs Mean: {} \n Min: {} \n Max: {} \n'.format(a.mean(),
# abs(a).mean(), a.min(), a.max()))


# '''
# Subtraction histogram plots
# '''
# for metal in [Os, U]:
#     '''
#     For each metal, saves histogram of the percent difference between the full-
#     and transmission-only measured intensities for each filter, as well as the
#     percent difference in full- and transmission-only subtractions between the
#     two filters. Also saves .png file of full subtraction
#     '''

#     plt.close()

#     if metal is Os:
#         nmet = 'Os'
#     elif metal is U:
#         nmet = 'U'

#     print('{}:'.format(nmet))

#     imAl, imAlT, imAlPh = forwardmodel(metal, Al)
#     imTi, imTiT, imTiPh = forwardmodel(metal, Ti)

#     print('Saving phase histogram...\n')
#     saveresults(imAlPh, '{}/{}_phase'.format(nmet, 'Al'))
#     saveresults(imTiPh, '{}/{}_phase'.format(nmet, 'Ti'))

#     print('imAlPh:')
#     print_stats(imAlPh)
#     print('imTiPh:')
#     print_stats(imTiPh)

#     # print('Saving histograms...')
#     # saveresults((imAl - imAlT) / imAl * 100,
#     #             '{}/{}_full_vs_T'.format(nmet, 'Al'))
#     # saveresults((imTi - imTiT) / imTi * 100,
#     #             '{}/{}_full_vs_T'.format(nmet, 'Ti'))

#     print('Energy subtraction histograms...')

#     sub = imAl - imTi
#     subT = imAlT - imTiT

#     E_effect = (sub - subT) / sub * 100

#     saveresults(E_effect, '{}/E_sub_full_vs_T'.format(nmet))
#     print_stats(E_effect)

#     # print('Saving image...\n')
#     # plt.imshow(sub, cmap=plt.cm.Greys_r)
#     # plt.axis('off')
#     # plt.tight_layout()
#     # plt.savefig('../Data/Spectra/Results/{}/sub_im.png'.format(nmet),
#     #             dpi=1000, bbox_inches='tight')
#     plt.close()
