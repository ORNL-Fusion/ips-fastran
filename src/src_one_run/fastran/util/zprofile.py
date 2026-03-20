from numpy import *
import scipy.optimize
from fastran.util.zinterp import zinterp

"""
  tanh0, ..
  Taken from Osborne pyped routine
"""
def tanh0(c,x,y,param=None):
    z = 2.*(c[0]-x)/c[1]
    pz = 1.
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  )
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    return yfit-y

def mtanh(c,x,y,param=None):
    z = 2.*(c[0]-x)/c[1]
    pz = 1. + c[4]*z
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  )
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    return yfit-y

def mtanh_parab(c,x,y,param=None):
    z = 2.*(c[0]-x)/c[1]
    pz = 1.
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  )
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    xped = c[0]-0.5*c[1]
    core = zeros(len(x))
    for k in range(len(x)):
        if x[k]<xped:
            core[k] = c[4]*(1.0-(x[k]/xped)**c[5])**c[6]
        else:
            core[k] = 0.
    yfit += core
    return yfit-y

def get_mtanh(c,x):
    z = 2.*(c[0]-x)/c[1]
    pz = 1. + c[4]*z
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  )
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    return yfit

def get_mtanh_parab(c,x):
    z = 2.*(c[0]-x)/c[1]
    pz = 1.
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  )
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    xped = c[0]-0.5*c[1]
    core = zeros(len(x))
    for k in range(len(x)):
        if x[k]<xped:
            core[k] = c[4]*(1.0-(x[k]/xped)**c[5])**c[6]
        else:
            core[k] = 0.
    yfit += core
    return yfit

class zprofile():
    func_maps = {"mtanh":mtanh, "tanh0":tanh0, "mtanh_parab":mtanh_parab}

    def __init__(self, x, y, func="mtanh", wped=0.1):
        self.x = x
        self.y = y
        self.func = self.func_maps[func]
        self.s = zinterp(self.x, self.y)
        self.wped = wped

    def info(self):
        pass

    def fit(self):
        wped_input = self.wped
        yped_input = self.s(1.0-self.wped)

        fit_coeff_input = [1.0-0.5*wped_input, wped_input, yped_input, 0.0, self.y[0]-yped_input, 1.5, 1.5]
        fit_coeff, success = scipy.optimize.leastsq(self.func, fit_coeff_input, args=(self.x, self.y))

        self.fit_coeff = fit_coeff

        self.xmid = fit_coeff[0]
        self.xwid = fit_coeff[1]
        self.yped = fit_coeff[2]
        self.yoff = fit_coeff[3]

    def __call__(self, xval):
#       rval = 0.0
#       rval = get_mtanh_parab(self.fit_coeff, xval)
        rval = get_mtanh(self.fit_coeff, xval)
        return rval

    def __getitem__(self, xvec):
        rvec = 0.0
        return rvec
