import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import laplace
from libtiff import TIFF

n = 512


def gauss2d(mu_x, mu_y, val=1, n=n):
    x = np.linspace(-10, 10, n)
    y = np.linspace(-10, 10, n)

    X, Y = np.meshgrid(x, y)
    sig = 1

    return val * np.exp(-((X + mu_x)**2 + (Y + mu_y)**2) / (2 * sig**2))


def circle(c_x=0, c_y=0, val=1, n=n):
    x = np.linspace(-10, 10, n)
    y = np.linspace(-10, 10, n)
    r = 1

    X, Y = np.meshgrid(x, y)
    Z = np.zeros(X.shape)
    Z[(Y - c_y)**2 <= r**2 - (X - c_x)**2] = val
    return Z


H2O = gauss2d(1, 6, 8) + \
    circle(-4.5, -4.5, 7) + \
    gauss2d(6, 0, 6) + \
    circle(-4.5, 4.5, 5) + \
    gauss2d(1, -6, 4) + \
    np.random.normal(0, 0.15, (n, n))


metal = circle(4.5, 4.5, 8) + \
    gauss2d(-1, -6, 7) + \
    circle(4.5, -4.5, 6) + \
    gauss2d(-6, 0, 5) + \
    gauss2d(-1, 6, 4) + \
    np.random.normal(0, 0.15, (n, n))

H2O -= H2O.min()
metal -= metal.min()

H2O /= H2O.mean()
metal /= metal.mean()

plt.figure(1)
plt.imshow(H2O)
plt.axis('off')
plt.savefig('../Data/Spectra/Results/new_phant/H2O.png',
            dpi=400, bbox_inches='tight')

plt.figure(2)
plt.imshow(metal)
plt.axis('off')
plt.savefig('../Data/Spectra/Results/new_phant/metal.png',
            dpi=400, bbox_inches='tight')


# np.save('../Data/Phantoms/h20.npy', H2O)
# np.save('../Data/Phantoms/metal.npy', metal)
