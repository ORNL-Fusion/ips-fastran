"""
 -----------------------------------------------------------------------
 model profiles
 -----------------------------------------------------------------------
"""
import numpy as np
import scipy.optimize
from fastran.util.zinterp import zinterp


def mtanh(c, x, param=None):
    """
    from Tom Osborne's pyped routines
    c0 = SYMMETRY POINT
    c1 = FULL WIDTH
    c2 = HEIGHT
    c3 = OFFSET
    c4 = SLOPE INNER
    c5 = QUADRADIC INNER
    c6 = CUBIC INNER
    """
    z  = 2. * ( c[0] - x ) / c[1]
    pz1  = 1.+ c[4] * z + c[5] * z * z + c[6] * z * z * z
    mth = 0.5 * ( pz1 + 1. ) * np.tanh(z) + 0.5 * (pz1 - 1.)
    y = 0.5*( (c[2] - c[3]) * mth + c[2] + c[3] )
    return y


def mtanh0(c, x, param=None):
    """
    from Tom Osborne's pyped routines
    c0 = SYMMETRY POINT
    c1 = FULL WIDTH
    c2 = HEIGHT
    c3 = OFFSET
    c4 = SLOPE INNER
    """
    z  = 2. * ( c[0] - x ) / c[1]
    pz1  = 1.+ c[4] * z
    mth = 0.5 * ( pz1 + 1. ) * np.tanh(z) + 0.5 * (pz1 - 1.)
    y = 0.5*( (c[2] - c[3]) * mth + c[2] + c[3] )
    return y


class profile_pedestal():
    def __init__(self, nx, xmid, xwidth, yped, ysep, yaxis, alpha, beta, ytop=0, ifit=0):
        x = np.arange(nx) / (nx - 1.)

        xped = xmid - 0.5 * xwidth
        xtop = xmid - xwidth

        a0 = (yped - ysep) / (np.tanh(2. * (1. - xmid) / xwidth) - np.tanh(2. * (xmid - 0.5 * xwidth - xmid) / xwidth))
        a1 = yaxis - ysep - a0 * (np.tanh(2. * (1. - xmid)/xwidth) - np.tanh(2. * (0. - xmid) / xwidth))
        if ytop > 0:
            yy = ysep + a0 * (np.tanh( 2. * (1. - xmid) / xwidth) - np.tanh(2. * (xtop - xmid) / xwidth))
            a1 = (ytop - yy) / (1. - (xtop / xped)**alpha )**beta
        y_edge = ysep + a0 * (np.tanh(2. * (1. - xmid) / xwidth) - np.tanh(2. * (x - xmid) / xwidth))

        y_core = np.zeros(nx)
        if yaxis > 0. or ytop > 0:
            for k, xval in enumerate(x):
                if xval < xped:
                    y_core[k] = a1 * (1. - (xval / xped)**alpha)**beta
                else:
                    y_core[k] = 0.

        y = y_edge + y_core

        self.nx = nx
        self.x = x
        self.y = y

        if ifit == 0:
            self.prof = zinterp(x, y)
        elif ifit == 1:
            def fit_func(c, x, y, param=None):
                return mtanh(c, x) - y
            c0 = [xmid, xwidth, yped, ysep, 0., 0., 0.] + 10 * [0.]
            self.c, success = scipy.optimize.leastsq(fit_func, c0, args=(x, y))
            self.prof = self.fit
        else:
            Exception('profile_pedestal: ifit option not valid')

    def __call__(self, xval):
        return self.prof(xval);

    def __getitem__(self, xvec):
        return 0.

    def fit(self, xval):
        return mtanh(self.c, xval)

    def info(self):
        pass


class profile_parab():
    def __init__(self, yaxis, ysep, alpha, beta):
        self.yaxis = yaxis
        self.ysep = ysep
        self.alpha = alpha
        self.beta = beta

    def __call__(self, xval):
        return (self.yaxis - self.ysep) * (1. - xval**self.alpha)**self.beta + self.ysep

    def __getitem__(self, xvec):
        return 0.

    def info(self):
        pass


class profile_spline():
    def __init__(self, x, y, yp):
        self.n = len(x)
        self.x = x
        self.y = y
        self.yp = yp

    def __call__(self, x):
        # from Numerical Recipes
        if x <= self.x[0]:
           return(self.y[0])
        elif x >= self.x[-1]:
           return(self.y[-1])
        else:
           klo = 0
           khi = self.n - 1
           while (khi - klo > 1):
              k=(khi + klo) >> 1
              if self.x[k] > x: khi = k
              else: klo = k
        h = self.x[khi] - self.x[klo]
        x1 = self.x[klo]
        x2 = self.x[khi]
        y1 = self.y[klo]
        y2 = self.y[khi]
        y1p = self.yp[klo]
        y2p = self.yp[khi]

        a = y1p * (x2 - x1) - (y2 - y1)
        b = -y2p * (x2 - x1) + (y2 - y1)
        t = (x - x1) / (x2 - x1)

        return (1. - t) * y1 + t * y2 + t * (1. - t) * (a * (1. - t) + b * t)


class profile_hat():
    def __init__(self, nx, r0, dr):
        x = np.arange(nx) / (nx - 1.)
        y = 1. - np.tanh(2. * (x - r0) / dr)

        self.nx = nx
        self.x = x
        self.y = y

        self.prof = zinterp(x, y)

    def info(self):
        pass

    def __call__(self,xval):
        return self.prof(xval)

    def __getitem__(self,xvec):
        return 0.


if __name__ == '__main__':
    pass
