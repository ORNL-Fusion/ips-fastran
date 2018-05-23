#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran instate file utility
 -----------------------------------------------------------------------
"""

from plasmastate import *
from numpy import *
from Namelist import Namelist
from zinterp import zinterp
from zmodelprof import profile_pedestal, profile_pedestal_smooth
import eped_io

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

def instate_model(fn_instate, instate_method):

    if instate_method == 'model':
        instate_model_default(fn_instate)
    else:
        eval("instate_model_%s(fn_instate)"%instate_method)

def instate_model_default(fn_instate):

    #-- read

    instate = Namelist(fn_instate)

    R0 = instate["instate"]["R0"][0]
    B0 = instate["instate"]["B0"][0]
    ip = instate["instate"]["ip"][0]
    nt = instate["instate"]["nbdry"][0]

    nrho = instate["instate"]["nrho"][0]
    xmid = instate["instate"]["xmid"][0]
    xwid = instate["instate"]["xwid"][0]

    ne0 = instate["instate"]["ne0"][0]
    neped  = instate["instate"]["neped"][0]
    nesep  = instate["instate"]["nesep"][0]
    alpha_ne = instate["instate"]["alpha_ne"][0]
    beta_ne = instate["instate"]["beta_ne"][0]

    te0 = instate["instate"]["te0"][0]
    teped = instate["instate"]["teped"][0]
    tesep = instate["instate"]["tesep"][0]
    alpha_te = instate["instate"]["alpha_te"][0]
    beta_te = instate["instate"]["beta_te"][0]

    ti0 = instate["instate"]["ti0"][0]
    tiped = instate["instate"]["tiped"][0]
    tisep = instate["instate"]["tisep"][0]
    alpha_ti = instate["instate"]["alpha_ti"][0]
    beta_ti = instate["instate"]["beta_ti"][0]

    zeff0 = instate["instate"]["zeff0"][0]

    omega0 = instate["instate"]["omega0"][0]
    alpha_omega = instate["instate"]["alpha_omega"][0]
    beta_omega  = instate["instate"]["beta_omega"][0]

    #-- pedestal

    if instate["instate"]["use_inped"][0] == 1:

        teped = instate["inped"]["teped_const"][0]
        teped*= ip**instate["inped"]["teped_ip"][0]
        teped*= B0**instate["inped"]["teped_bt"][0]
        teped*= neped**instate["inped"]["teped_neped"][0]

        wped = instate["inped"]["wped_const"][0]
        wped*= ip**instate["inped"]["wped_ip"][0]
        wped*= B0**instate["inped"]["wped_bt"][0]
        wped*= neped**instate["inped"]["wped_neped"][0]

        xmid = 1.0-0.5*wped
        xwid = wped

        niped, nzped = eped_io.get_ni(neped,zeff=zeff0)
        tiped = teped * neped / (niped+nzped)
        ti0 = te0*tiped/teped

        nesep = neped * instate["inped"]["nesep"][0]

        tesep_model =  instate["inped"]["tesep"][0]
        if tesep_model > 0:
            tesep = tesep_model
        else:
            tesep = abs(tesep_model)*teped

        tisep_model =  instate["inped"]["tisep"][0]
        if tisep_model > 0:
            tisep = tisep_model
        else:
            tisep = abs(tisep_model)*tiped

        print 'teped, tiped = ', teped, tiped
        print 'niped, nzped = ', niped, nzped
        print 'nesep = ', nesep
        print 'tesep, tisep = ', tesep, tisep

    #-- make profiles

    rho = linspace(0.0,1.0,nrho)
    ne = profile_pedestal(nrho, xmid, xwid, neped, nesep, ne0, alpha_ne, beta_ne)(rho)
    te = profile_pedestal(nrho, xmid, xwid, teped, tesep, te0, alpha_te, beta_te)(rho)
    ti = profile_pedestal(nrho, xmid, xwid, tiped, tisep, ti0, alpha_ti, beta_ti)(rho)

    zeff = nrho*[zeff0]
    omega = omega0*(1.-rho**alpha_omega)**beta_omega
    j_tot = 1.0-rho**2 + 0.05

    #-- boundary

    if instate["instate"]["model_shape"][0] == 1:
        rb, zb, rlim, zlim = \
            set_shape(
                R0 = instate["instate"]["r0"][0],
                a0 = instate["instate"]["a0"][0],
                kappa = instate["instate"]["kappa"][0],
                delta = instate["instate"]["delta"][0],
                nt = instate["instate"]["nbdry"][0])
    else:
        rb, zb, rlim, zlim = \
            instate["instate"]["rbdry"], \
            instate["instate"]["zbdry"], \
            instate["instate"]["rlim"], \
            instate["instate"]["zlim"]

    #-- clean up

    del instate["instate"]["nrho"]
    del instate["instate"]["xmid"]
    del instate["instate"]["xwid"]
    del instate["instate"]["ne0"]
    del instate["instate"]["neped"]
    del instate["instate"]["nesep"]
    del instate["instate"]["alpha_ne"]
    del instate["instate"]["beta_ne"]
    del instate["instate"]["te0"]
    del instate["instate"]["teped"]
    del instate["instate"]["tesep"]
    del instate["instate"]["alpha_te"]
    del instate["instate"]["beta_te"]
    del instate["instate"]["ti0"]
    del instate["instate"]["tiped"]
    del instate["instate"]["tisep"]
    del instate["instate"]["alpha_ti"]
    del instate["instate"]["beta_ti"]
    del instate["instate"]["omega0"]
    del instate["instate"]["alpha_omega"]
    del instate["instate"]["beta_omega"]
    del instate["instate"]["zeff0"]

    #-- put profiles to instate

    instate["instate"]["nrho"] = [nrho]
    instate["instate"]["rho"] = rho

    instate["instate"]["ne"] = ne
    instate["instate"]["te"] = te
    instate["instate"]["ti"] = ti
    instate["instate"]["omega"] = omega
    instate["instate"]["zeff"] = zeff

    instate["instate"]["density_alpha"] = zeros(nrho)
    instate["instate"]["walpha"] = zeros(nrho)
    instate["instate"]["density_beam"] = zeros(nrho)
    instate["instate"]["wbeam"] = zeros(nrho)

    instate["instate"]["j_tot"] = j_tot

    instate["instate"]["nbdry"] = [len(rb)]
    instate["instate"]["rbdry"] = rb
    instate["instate"]["zbdry"] = zb
    instate["instate"]["nlim"] = [len(rlim)]
    instate["instate"]["rlim"] = rlim
    instate["instate"]["zlim"] = zlim

    instate.write(fn_instate)
