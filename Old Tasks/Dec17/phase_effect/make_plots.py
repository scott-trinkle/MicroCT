'''
SPECTRUM PLOT FOR GROUP MEETING
'''
# plt.plot(Al.E, Al.I0, label='1 mm Al, 3 mrad')
# plt.plot(Ti.E, Ti.I0, label='1 mm Ti, 2 mrad')
# plt.xlabel('E [keV]')
# plt.ylabel('Intensity [photons/s/0.1% BW/$mrad^2$')
# plt.title(r'$I_0^{(j)}$')
# plt.legend()
# plt.savefig('../Notes/group_meeting_11_15/figs/spectra.png', dpi=400)

'''
Laplacian of Phi - histogram
THIS GOES INSIDE "forwardmodel" in phasemagnitude.py
'''
# Save laplacian of phi histograms
# plt.close()
# plt.hist(lap_phi_E.flatten(), bins=500)
# plt.xlabel(r'$\nabla^2 \phi(E)$')
# plt.title('Spectrum: {}\n Metal: {}'.format(nspect, nmet))
# plt.savefig(
#     '../Data/Spectra/data_plots/lap_{}_{}.png'.format(nspect, nmet), dpi=400)

'''
LAPLACIAN of n_{a,i} PLOT
'''

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


'''
F1 AND F2 PLOT
'''

# plt.close()
# axf1 = plt.subplot(2, 1, 1)
# axf2 = plt.subplot(2, 1, 2)
# labels = ['H2O', 'Os', 'U']

# for i, mat in enumerate([H2O, Os, U]):
#     axf1.plot(mat.E, mat.f1, label=labels[i])
#     axf1.legend()
#     axf2.plot(mat.E, mat.f2, label=labels[i])
#     axf2.legend()

# axf1.set_title('f1 vs. E')
# axf1.set_xlabel('E [keV]')
# axf1.set_ylabel('f1')
# axf2.set_title('f2 vs. E')
# axf2.set_xlabel('E [keV]')
# axf2.set_ylabel('f2')
# plt.tight_layout()

# plt.savefig('../Data/Spectra/f1_f2.png', dpi=400, bbox_inches='tight')
