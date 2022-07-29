# Luce shape
from omfit_classes.fluxSurface import boundaryShape
#from omfit.classes.fluxSurface import boundaryShape

# simple analytic
def set_shape(R0, a0, kappa, delta, nt, dlim = 0.05):
    t = linspace(0.0, 2.*pi, nt)
    rb = R0 + a0*cos(t+delta*sin(t))
    zb = kappa*a0*sin(t)

    rmax = max(rb) + dlim
    rmin = min(rb) - dlim
    zmax = max(zb) + dlim
    zmin = min(zb) - dlim
    rlim = [ rmax, rmin, rmin, rmax, rmax ]
    zlim = [ zmax, zmax, zmin, zmin, zmax ]

    return rb, zb, rlim, zlim


