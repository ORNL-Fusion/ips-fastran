from numpy import *

kB = 1.6022e-19
e  = 1.6022e-19
mp = 1.6726e-27
me = 9.1094e-31

def rhostar(A,Ti,B,a):
# keV, Tesla, m
#   return (2.*(A*mp)*kB*Ti*1.0e3)**0.5/(e*B)/a
    return 0.00456933164275*(A*Ti)**0.5/B/a

def nustar_(R,r,ne,te,zeff,q):
# m, 10^19/m^3, keV

    neCGS = ne*1.e13
    teEV = te*1.e3
    lnlamda = 24.0-0.5*log(neCGS)+log(teEV);
    taue = 3.44e5*teEV**1.5/(neCGS*lnlamda*zeff)
    xnue = 0.75*pi**0.5/taue
    vthe = (kB*teEV/me)**0.5
    eps = r/R
    return xnue/vthe*eps**-1.5*q*R, xnue #lnlamda #vthe*sqrt(2)

def nustar(R,r,ne,te,zeff,q):
# m, 10^19/m^3, keV

    vthm1=5.9269e7
    tauenorm=3.443e-8

    neCGS = ne
    teEV = te*1.e3
    lnlamda = 16.7-0.5*log(neCGS/(0.001*teEV)**2)
    taue = tauenorm*teEV**1.5/(neCGS*zeff*lnlamda)
    xnue = 1./taue
    vthe = vthm1*teEV**0.5
    eps = r/R
    return sqrt(2)*xnue/vthe*eps**-1.5*q*R*1.0e2,xnue #lnlamda #vthe
