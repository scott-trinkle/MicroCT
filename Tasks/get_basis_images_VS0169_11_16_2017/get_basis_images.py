import numpy as np
import matplotlib.pyplot as plt
from libtiff import TIFF
from scipy.interpolate import interp1d

# Uranium
# Nominal density: 1.892e1 g/cm3
# E keV             u/p  cm2/g

U = np.array([[1.708047E+01,  4.4787E+01],
              [1.714913E+01,  4.4345E+01],
              [1.725213E+01,  1.0327E+02],
              [1.750963E+01,  9.9419E+01],
              [1.759961E+01,  9.8122E+01]])

# Osmium
# Nominal density: 2.253e1 g/cm3
# E kev           u/p cm2/g
Os = np.array([[1.065348E+01,  8.5388E+01],
               [1.081655E+01, 8.2028E+01],
               [1.086003E+01,  8.1166E+01],
               [1.092525E+01,  2.0788E+02],
               [1.103212E+01,  2.0250E+02]])

# Water
H2O = np.array([[7.902609E+00,  1.027736E+01],
                [8.447890E+00,  8.464951E+00],
                [9.030794E+00,  6.971115E+00],
                [9.653919E+00,  5.721777E+00],
                [1.032004E+01,  4.696412E+00],
                [1.103212E+01,  3.855846E+00],
                [1.179334E+01,  3.174232E+00],
                [1.260708E+01,  2.621632E+00]])

U_interp_func = interp1d(U[:, 0], U[:, 1], kind='linear')
Os_interp_func = interp1d(Os[:, 0], Os[:, 1], kind='linear')
H2O_interp_func = interp1d(H2O[:, 0], H2O[:, 1], kind='linear')

U_1712 = U_interp_func(17.12)
U_1722 = U_interp_func(17.22)
H2O_1712 = H2O_interp_func(17.12)
H2O_1722 = H2O_interp_func(17.22)

Os_1082 = Os_interp_func(10.82)
Os_1092 = Os_interp_func(10.92)
H2O_1082 = H2O_interp_func(10.82)
H2O_1092 = H2O_interp_func(10.92)
