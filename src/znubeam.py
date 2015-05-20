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

########################################################################
# io

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

    f_ion = ns_ion/sum(ns_ion,axis=0)
    f_imp = ns_imp/sum(ns_imp,axis=0)

    #print f_ion
    #print f_imp

    # DLG : Get beam density on rho grid
    # This is a terrible terrible hack.
    # Assumes dn/drho=0 at boundaries.
    # Should use the plasma-state API to get
    # the nbeami data at the rho grid desired.
    # - OR - 
    # looks like JM has some cell2node functions
    # that could do part of this.

    _rho_nbi = ps["rho_nbi"][:]
    _drho_nbi = _rho_nbi[1]-_rho_nbi[0]
    _rho_nbi_dm1 = _rho_nbi[0:nrho_nbi-1]+_drho_nbi/2
    _rho_nbi_full = concatenate([[_rho_nbi[0]-_drho_nbi],_rho_nbi_dm1,[_rho_nbi[-1]+_drho_nbi]])
    _density_beam_full = concatenate([[density_beam[0]],density_beam,[density_beam[-1]]])
    _density_beam = interp1d(_rho_nbi_full,_density_beam_full)

    for i in range(nrho-1):

        # print i, f_ion[:,i], f_imp[:,i]
        
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

        #zne_adj = zne_adj - 1.0*density_beam[i]
        #zzne_adj = zzne_adj - 1.0**2*density_beam[i]

        # DLG : This is part of the terrible hack
        zne_adj = zne_adj - 1.0*_density_beam(rho[i])
        zzne_adj = zzne_adj - 1.0**2*_density_beam(rho[i])

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
