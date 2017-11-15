import numpy as np
import matplotlib.pyplot as plt
from phasefunctions import Material, Spectrum, laplace_spectrum, match_arrays, calc_lap_phi, calc_T
from libtiff import TIFF
from scipy.ndimage.filters import laplace
from scipy.ndimage import zoom


# Reading in energy spectra
Al = Spectrum('../Data/Spectra/Al_1mm_3mrad.txt')
Ti = Spectrum('../Data/Spectra/Ti_1mm_2mrad.txt')

# Reading in f1 and f2 values
Os = Material('../Data/Spectra/f1_Os.txt', '../Data/Spectra/f2_Os.txt')
U = Material('../Data/Spectra/f1_U.txt', '../Data/Spectra/f2_U.txt')
H2O = Material('../Data/Spectra/f1_H2O_interp.txt',
               '../Data/Spectra/f2_H20.txt')

# SPECTRUM PLOT FOR GROUP MEETING
# plt.plot(Al.E, Al.I0, label='1 mm Al, 3 mrad')
# plt.plot(Ti.E, Ti.I0, label='1 mm Ti, 2 mrad')
# plt.xlabel('E [keV]')
# plt.ylabel('Intensity [photons/s/0.1% BW/$mrad^2$')
# plt.title(r'$I_0^{(j)}$')
# plt.legend()
# plt.savefig('../Notes/group_meeting_11_15/figs/spectra.png', dpi=400)


# Constants
r_e = 2.818e-15  # m [Jacobsen, Kirz, Howells chapter]
R2 = 0.30  # m [.setup files]

Na = 6.02214e23  # atoms/mole, https://en.wikipedia.org/wiki/Avogadro_constant

Os.A = 190.23  # g / mol    https://en.wikipedia.org/wiki/Osmium
U.A = 238.02891  # g / mol  https://en.wikipedia.org/wiki/Uranium
H2O.A = 18.01488  # g / mol https://en.wikipedia.org/wiki/Molar_mass

pix_size = 1.24e-6  # m [.setup files]
width = 1920 * pix_size  # m

# Calculating number densities with formular from [Jacobsen, Kirz, Howells]
# Densities from the above cited wiki pages.
# Note these are actually the number density LINE INTEGRALS in
# particles / cm^2.
H2O.ndens = 1 * Na / H2O.A * width * 100e3  # particles / m^2
Os.ndens = 22.59 * Na / Os.A * width * 100e3  # particles / m^2
U.ndens = 19.1 * Na / U.A * width * 100e3  # particles / m^2

# Reading in image files to calculate "Laplacian Spectrums"
Al_phant = TIFF.open('../Data/Projection_Tifs/Al_1080.tif', 'r').read_image()
No_phant = TIFF.open('../Data/Projection_Tifs/No_1080.tif', 'r').read_image()

# Al_phant = zoom(Al_phant, (1 / 10, 1 / 10), order=0)
# No_phant = zoom(No_phant, (1 / 10, 1 / 10), order=0)

# Al_phant = np.load('../Data/Phantoms/h20.npy')
# No_phant = np.load('../Data/Phantoms/metal.npy')

H2O.phant = Al_phant / Al_phant.mean() * H2O.ndens
Os.phant = No_phant / No_phant.mean() * Os.ndens
U.phant = No_phant / No_phant.mean() * U.ndens

print('\nCalculating laplacians:')
print('H2O...')
H2O.lap = laplace(H2O.phant)
print('Os...')
Os.lap = laplace(Os.phant)
print('U...\n')
U.lap = laplace(U.phant)


def forwardmodel(metal=Os, spect=Al, H2O=H2O):
    '''
    Matches all array sizes, calculates laplacian of phi and transmission factor,
    Integrates over energy to calculate the full measured intensity, Im_full
    and the transmission-only intensity Im_trans
    '''

    if metal is Os:
        nmet = 'Os'
    elif metal is U:
        nmet = 'U'
    if spect is Al:
        nspect = 'Al'
    elif spect is Ti:
        nspect = 'Ti'

    E, H2O, spect, metal = match_arrays(H2O, spect, metal)

    print('Calculating lap_phi_E...')
    lap_phi_E = calc_lap_phi(spect, H2O, metal)

    plt.close()
    plt.hist(lap_phi_E.flatten(), bins=500)
    plt.xlabel(r'$\nabla^2 \phi(E)$')
    plt.title('Spectrum: {}\n Metal: {}'.format(nspect, nmet))
    plt.savefig(
        '../Data/Spectra/data_plots/lap_{}_{}.png'.format(nspect, nmet), dpi=400)

#     print('Calculating T_E...')
#     T_E = calc_T(spect, H2O, metal)

#     print('Calculating Im...')
#     Im_full = np.array([E[i] * spect.I0_int[i] *
#                         T_E[:, :, i] * (1 + R2 / spect.K_int[i] *
#                                         lap_phi_E[:, :, i]) for i in range(E.size)]).sum(axis=0)
#     Im_trans = np.array([E[i] * spect.I0_int[i] * T_E[:, :, i]
#                          for i in range(E.size)]).sum(axis=0)

#     # Im_phase = np.array([E[i] * spect.I0_int[i] *
#     #                      (1 + R2 / spect.K_int[i] *
#     # lap_phi_E[:, :, i]) for i in range(E.size)]).sum(axis=0)

#     # ph_factor = np.array([(1 + R2 / spect.K_int[i] * lap_phi_E[:, :, i])
#     #                       for i in range(E.size)])

#     return Im_full, Im_trans


# def saveresults(result, fn):
#     '''
#     Saves a histogram of the results
#     '''
#     plt.close()
#     plt.hist(result.flatten(), bins=500)
#     if fn is not None:
#         plt.title(fn)
#         plt.xlabel('Percent Difference')
#         plt.savefig('../Data/Spectra/Results/new_phant/' +
#                     fn + '.png', dpi=400)
#         plt.close()
#     else:
#         plt.show()

#                                 LAPLACIAN of n_{a,i} PLOT


# fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

# axes[0].hist(H2O.lap.flatten(), bins=500, label='H2O')
# axes[0].legend()


# axes[1].hist(Os.lap.flatten(), bins=500, label='Os')
# axes[1].legend()


# axes[2].hist(U.lap.flatten(), bins=500, label='U')
# axes[2].legend()
# axes[2].set_xlabel('[m$^{-4}$]')

# fig.suptitle(r'$\nabla^2 \int n_{a,i}(\vec{x}) dl$')

# plt.savefig('../Data/Spectra/data_plots/laplacians.png', dpi=400)


# #                               F1 AND F2 PLOT

# # plt.close()
# # axf1 = plt.subplot(2, 1, 1)
# # axf2 = plt.subplot(2, 1, 2)
# # labels = ['H2O', 'Os', 'U']

# # for i, mat in enumerate([H2O, Os, U]):
# #     axf1.plot(mat.E, mat.f1, label=labels[i])
# #     axf1.legend()
# #     axf2.plot(mat.E, mat.f2, label=labels[i])
# #     axf2.legend()

# # axf1.set_title('f1 vs. E')
# # axf1.set_xlabel('E [keV]')
# # axf1.set_ylabel('f1')
# # axf2.set_title('f2 vs. E')
# # axf2.set_xlabel('E [keV]')
# # axf2.set_ylabel('f2')
# # plt.tight_layout()

# # plt.savefig('../Data/Spectra/f1_f2.png', dpi=400, bbox_inches='tight')


# #                        Subtraction histogram plots
# for metal in [Os, U]:
#     '''
#     For each metal, saves histogram of the percent difference between the full-
#     and transmission-only measured intensities for each filter, as well as the
#     percent difference in full- and transmission-only subtractions between the
#     two filters. Also saves .png file of full subtraction
#     '''

#     if metal is Os:
#         nmet = 'Os'
#     elif metal is U:
#         nmet = 'U'

#     print('{}:'.format(nmet))

#     imAl, imAlT, imAlPh = forwardmodel(metal, Al)
#     imTi, imTiT, imTiPh = forwardmodel(metal, Ti)

#     print('Saving phase histogram...')
#     saveresults((imAlPh - imTiPh) / imAlPh * 100, '{}/{}_full_vs_ph')

#     print('Saving histograms...')
#     saveresults((imAl - imAlT) / imAl * 100,
#                 '{}/{}_full_vs_T'.format(nmet, 'Al'))
#     saveresults((imTi - imTiT) / imTi * 100,
#                 '{}/{}_full_vs_T'.format(nmet, 'Ti'))

#     sub = imAl - imTi
#     subT = imAlT - imTiT

#     saveresults((sub - subT) / sub * 100, '{}/E_sub_full_vs_T'.format(nmet))

#     print('Saving image...\n')
#     plt.imshow(sub, cmap=plt.cm.Greys_r)
#     plt.axis('off')
#     plt.tight_layout()
#     plt.savefig('../Data/Spectra/Results/{}/sub_im.png'.format(nmet),
#                 dpi=1000, bbox_inches='tight')
#     plt.close()
