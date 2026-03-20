import numpy as np
# Luce shape
from omfit_classes.fluxSurface import boundaryShape


# simple analytic
def set_shape(R0, a0, kappa, delta, nt):
    t = np.linspace(0., 2. * np.pi, nt)
    rb = R0 + a0 * np.cos(t + delta * np.sin(t))
    zb = kappa * a0 * np.sin(t)

    return rb, zb

# box limiter
def set_limiter(rb, zb, dlim=0.05):
    rmax = np.max(rb) + dlim
    rmin = np.min(rb) - dlim
    zmax = np.max(zb) + dlim
    zmin = np.min(zb) - dlim
    rlim = [rmax, rmin, rmin, rmax, rmax]
    zlim = [zmax, zmax, zmin, zmin, zmax]

    return rlim, zlim
