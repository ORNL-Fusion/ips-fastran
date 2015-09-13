#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 utils for nubeam IO
 JM
 -----------------------------------------------------------------------
"""

import sys,os,glob,pickle,shutil
from numpy import * 
from scipy.interpolate import interp1d

import Namelist
from zplasmastate import plasma_state_file
from zefitutil import readg
from plasmastate import *

########################################################################
# io

def ps_interp1d(rho_, f_, rho):

    # -----------------------------------------
    # Interpolates from the dm1 plasma-state grid
    # Underscore vars are on the dm1 (minus 1) plasma-state grid
    # This assumes df/dRho=0 at boundaries

    drho = rho_[1]-rho_[0]
    rhoCentered = rho_[0:-1]+drho/2
    rhoExtend = concatenate([[rhoCentered[0]-drho],rhoCentered,[rhoCentered[-1]+drho]])
    fExtend = concatenate([[f_[0]],f_,[f_[-1]]])
    f_interp = interp1d(rhoExtend,fExtend)

    return f_interp(rho)


def update_ps_profile(f_state,f_eqdsk):

    #------------------------------------------------------------------- 
    # read geqdsk and plasma state file

    geq = readg(f_eqdsk) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']
    ps  = plasma_state_file(f_state,r0=r0,b0=b0,ip=ip)

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27

    z_spec = [round(x) for x in ps["qatom_S"][1:]/ps_xe ]
    a_spec = [round(x) for x in ps["m_S"][1:]/ps_mp ]
    n_spec = len(z_spec)

    n_imp = len(ps.data.dimensions["dim_nspec_imp0"])
    n_ion = n_spec-n_imp

    z_ion = z_spec[0:n_ion]
    a_ion = a_spec[0:n_ion]
    z_imp = z_spec[n_ion:n_imp+n_ion]
    a_imp = a_spec[n_ion:n_imp+n_ion]

    ns_ion = ps["ns"][1:n_ion+1]
    ns_imp = ps["ns"][n_ion+1:n_imp+n_ion+1]

    nrho  = ps.nrho
    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]
    nrho_nbi = len(ps.data.dimensions["dim_nrho_nbi"])
    density_beam = ps["nbeami"][0,:]
    zeff  = ps["Zeff"][:]

    #------------------------------------------------------------------
    # density

    density_ion = {}
    for k in range(n_ion):
        density_ion[k] = zeros(nrho-1)

    density_imp = {}
    for k in range(n_imp):
        density_imp[k] = zeros(nrho-1)

    density_th = zeros(nrho)

    # Hack to avoid zero density problems

    min_density = 1e-5
    ns_imp_zero_ii = ns_imp < min_density
    ns_imp[ns_imp_zero_ii] = min_density

    f_ion = ns_ion/sum(ns_ion,axis=0)
    f_imp = ns_imp/sum(ns_imp,axis=0)

    # The plasma-state API inter1d function does not work.
    #_ps = PlasmaState("ips",1)
    #_ps.read(f_state)
    #density_beam_at_rho = _ps.interp1d("nbeami",0,rho) 

    density_beam_at_rho = ps_interp1d(ps["rho_nbi"][:],ps["nbeami"][0,:],rho)

    for i in range(nrho-1):

        a=0; b=0; c=0; d=0
        for k in range(n_imp):
            b = b+f_imp[k,i]*z_imp[k]
            d = d+f_imp[k,i]*z_imp[k]*z_imp[k]
        for k in range(n_ion):
            a = a+f_ion[k,i]*z_ion[k]
            c = c+f_ion[k,i]*z_ion[k]*z_ion[k]

        zne_adj = ne[i]
        zzne_adj = ne[i]*zeff[i]

        # depletion due to beam ions

        zne_adj = zne_adj - 1.0*density_beam_at_rho[i]
        zzne_adj = zzne_adj - 1.0**2*density_beam_at_rho[i]

        # effective main ion and impurity densities

        nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
        nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

        for k in range(n_ion):
            density_ion[k][i] = f_ion[k,i]*nion
        for k in range(n_imp):
            density_imp[k][i] = f_imp[k,i]*nimp

        for k in range(n_ion):
            density_th[i] = density_th[i] + density_ion[k][i]
        for k in range(n_imp):
            density_th[i] = density_th[i] + density_imp[k][i]

    for k in range(n_ion):
        #print ps["ns"][k+1][:]
        ps["ns"][k+1] = density_ion[k]
        #print ps["ns"][k+1][:]

    ps.close()
