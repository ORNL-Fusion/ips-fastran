#!/usr/bin/env python

"""
 ------------------------------------------------------------------------------
  metric from efit
 ------------------------------------------------------------------------------
"""

import sys,os,re,pickle,glob
from numpy import *
from scipy.interpolate import interp1d
import Namelist

from zinterp import zinterp

def zinmetric(ps,r0,b0,ip,nrho=51,iwrt=False):

    psi     = ps["psipol"][:]/ps["psipol"][-1]
    rho     = sqrt(ps["phit"][:]/ps["phit"][-1])
    nrho    = len(rho)
    rhob    = (ps["phit"][-1]/pi/abs(b0))**0.5
    rhopsi  = zinterp(psi,rho)
    ipol    = ps["g_eq"][:]/(r0*b0)
    ipol    = abs(ipol)
    volp    = 4.0*pi*pi*rho*rhob*r0/ipol/(r0*r0*ps["gr2i"][:])
    g11     = volp*ps["grho2"][:]*rhob**2
    g22     = r0*volp/(4.0*pi*pi)*ps["grho2r2i"][:]*rhob**2
    g33     = r0*r0*ps["gr2i"][:]
    gradrho = ps["grho1"][:]*rhob
    area    = ps["surf"][:]
    rmajor  = ps["Rmajor_mean"][:]
    rminor  = ps["rMinor_mean"][:]
    shift   = rmajor-r0
    kappa   = ps["elong"][:]
    delta   = ps["triang"][:]
    pmhd    = ps["P_eq"][:]
    qmhd    = ps["q_eq"][:]

    nc1 = ps["gncfh"][:]
    gb1 = ps["gb1"][:]
    gb2 = ps["gb2"][:]
    Bmax = ps["B_surfMax"][:]
    hfac1 = gb1/Bmax
    hfac2 = gb2/Bmax**2

    bp2 = ps['gb2'][:] - ps['g_eq'][:]**2 * ps['gr2i'][:]

    metric = Namelist.Namelist()

    metric['inmetric']['Ip'     ] = [ip]
    metric['inmetric']['bcentr' ] = [b0]
    metric['inmetric']['rmajor' ] = [r0]
    metric['inmetric']['aminor' ] = [rminor[-1]]
    metric['inmetric']['kappa'  ] = [kappa [-1]]
    metric['inmetric']['delta'  ] = [delta [-1]]

    metric['inmetric']['nrho'   ] = [nrho]
    metric['inmetric']['rhob'   ] = [rhob]
    metric['inmetric']['rho'    ] = rho
    metric['inmetric']['volp'   ] = volp
    metric['inmetric']['ipol'   ] = ipol
    metric['inmetric']['g11'    ] = g11
    metric['inmetric']['g22'    ] = g22
    metric['inmetric']['g33'    ] = g33
    metric['inmetric']['gradrho'] = gradrho
    metric['inmetric']['area'   ] = area
    metric['inmetric']['rminor' ] = rminor
    metric['inmetric']['rmajor' ] = rmajor
    metric['inmetric']['shift'  ] = shift
    metric['inmetric']['kappa'  ] = kappa
    metric['inmetric']['delta'  ] = delta
    metric['inmetric']['pmhd'   ] = pmhd
    metric['inmetric']['qmhd'   ] = qmhd
    metric['inmetric']['er'     ] = zeros(nrho) #<======
    metric['inmetric']['nc1'    ] = nc1
    metric['inmetric']['hfac1'  ] = hfac1
    metric['inmetric']['hfac2'  ] = hfac2
    metric['inmetric']['psi'    ] = ps["psipol"][:]
    metric['inmetric']["vol"    ] = ps["vol"][:]
    metric['inmetric']["gr2i"   ] = ps["gr2i"][:]
    metric['inmetric']["bp2"    ] = bp2

    return metric
