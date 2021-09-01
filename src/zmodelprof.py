"""
 -----------------------------------------------------------------------
 model profile with tanh pedestal
 -----------------------------------------------------------------------
"""

from numpy import *
import scipy.optimize
import zinterp

class profile_pedestal():
    def __init__(self, nx, xmid, xwidth, yped, ysep, yaxis, alpha, beta, ytop=0, ifit=0):
        xped = xmid-xwidth/2
        xtop = xmid-xwidth

        a0 = (yped-ysep)/(tanh(2.0*(1-xmid)/xwidth)-tanh(2.0*(xmid-0.5*xwidth-xmid)/xwidth))
        a1 = yaxis - ysep - a0*(tanh(2.0*(1-xmid)/xwidth)-tanh(2.0*(0.0-xmid)/xwidth))
        if ytop > 0.0:
            yy = ysep + a0*(tanh(2.0*(1-xmid)/xwidth)-tanh(2.0*(xtop-xmid)/xwidth))
            a1 = (ytop-yy)/(1.0-(xtop/xped)**alpha)**beta

        x = arange(nx)/(nx-1.0)
        y_edge = ysep + a0*(tanh(2.0*(1-xmid)/xwidth)-tanh(2.0*(x-xmid)/xwidth))

        y_core = zeros(nx)
        if yaxis > 0.0 or ytop > 0:
            for k,xval in enumerate(x):
                if xval < xped: y_core[k] = a1*(1.0-(xval/xped)**alpha)**beta
                else: y_core[k] = 0.0

        y = y_edge+y_core

        self.nx = nx
        self.x  = x
        self.y  = y

        if ifit == 0:
            self.prof = zinterp.zinterp(x, y)
        elif ifit == 1:
            c0 = [xmid, xwidth, yped, ysep, 0.0, 0.0, 0.0] + 10*[0.0]
            self.c, success = scipy.optimize.leastsq(self.fit_func, c0, args=(x, y, self.tanhg))
            self.prof = self.call_tanhg
        else:
            Exception ('profile_pedestal: ifit option not valid')

    def __call__(self, xval):
        return self.prof(xval);

    def __getitem__(self, xvec):
        return 0.0

    def info(self):
        pass

    def tanhg(self, c, x, param=None):
        """
          Taken from Osborne pyped routine
          c[0] = SYMMETRY POINT
          c[1] = FULL WIDTH
          c[2] = HEIGHT
          c[3] = OFFSET
        """
        z = 2.*(c[0]-x)/c[1]
        pz1  = 1.+ c[4]*z + c[5]*z*z + c[6]*z*z*z
        pz2  = 1.+ c[7]*z
        mth = 0.5*( ( pz1 + pz2 )*tanh(z) + pz1 - pz2  )
        y = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
        return y

    def fit_func(self, c, x, y, func, param=None):
        rval = func(c,x)-y
        return rval

    def call_tanhg(self, xval):
        return self.tanhg(self.c, xval)

class profile_spline():
    def __init__(self, x, y, yp):
        self.n = len(x)
        self.x = x
        self.y = y
        self.yp = yp

    def __call__(self, x):
        if x <= self.x[0]:
           return(self.y[0])
        elif x >= self.x[-1]:
           return(self.y[-1])
        else:
           klo=0
           khi=self.n-1
           while (khi-klo > 1):
              k=(khi+klo) >> 1
              if self.x[k]>x: khi=k
              else: klo=k
        h = self.x[khi]-self.x[klo]
        x1 = self.x[klo]
        x2 = self.x[khi]
        y1 = self.y[klo]
        y2 = self.y[khi]
        y1p = self.yp[klo]
        y2p = self.yp[khi]

        a = y1p*(x2-x1)-(y2-y1)
        b =-y2p*(x2-x1)+(y2-y1)
        t = (x-x1)/(x2-x1)

        return (1-t)*y1 + t*y2 + t*(1.0-t)*(a*(1.0-t)+b*t)

class profile_hat():
    def __init__(self, nx, r0, dr):
        x = arange(nx)/(nx-1.0)
        y = 1-tanh( 2.0*(x-r0)/dr )

        self.nx = nx
        self.x  = x
        self.y  = y

        self.prof = zinterp.zinterp(x,y)

    def info(self):
        pass

    def __call__(self,xval):
        return self.prof(xval);

    def __getitem__(self,xvec):
        return 0.0

if __name__ == "__main__":
    x0 = [0.0,0.6,0.8]
    y0 = [0.3,1.0,0.2]
    yp = [1.0,0.0,-1.0]
    spl = profile_spline(x=x0,y=y0,yp=yp)

    x = linspace(0.0,0.8,num=100)
    y = [ spl(xx) for xx in x ]
    print (y)
