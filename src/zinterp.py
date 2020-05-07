
"""
 -----------------------------------------------------------------------
 interpolation/fit wrapper
 -----------------------------------------------------------------------
"""

from numpy import *
from scipy.interpolate import interp1d, splrep, splev

class zinterp():

    def __init__(self, x, y, s=0):
        self.x = x
        self.y = y
        self.dat = splrep(x, y, s=s)

    def info(self):
        pass

    def __call__(self, xval, der=0):
        rval = splev(xval, self.dat, der=der)
        return rval

    def __getitem__(self, xvec, der=0):
        rvec = splev(xvec, self.dat, der=der)
        return rvec

if __name__ == "__main__":
    pass
