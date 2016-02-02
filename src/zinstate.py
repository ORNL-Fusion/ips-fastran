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

#if instate["nrho"][0] != ps["nrho"]

def instate2ps(fn_instate,ps):

    #------------------------------------------------------------------
    # from instate

    instate = Namelist(fn_instate)["instate"]

    for key in instate.keys(): instate[key] = array(instate[key])

    nrho  = instate["nrho"][0]
    rho   = instate["rho"]

    n_ion = instate["n_ion"][0]
    z_ion = instate["z_ion"]
    a_ion = instate["a_ion"]
    f_ion = instate["f_ion"]

    n_imp = instate["n_imp"][0]
    z_imp = instate["z_imp"]
    a_imp = instate["a_imp"]
    f_imp = instate["f_imp"]

    ne = instate["ne"]
    te = instate["te"]
    ti = instate["ti"]
    omega = instate["omega"]
    zeff = instate["zeff"]

    ne[0] = ne[1]
    te[0] = te[1]
    ti[0] = ti[1]
    omega[0] = omega[1]
    zeff[0] = zeff[1]

    r0 = instate["r0"][0]
    b0 = instate["b0"][0]

    #------------------------------------------------------------------
    # put zeros if not defined

    for key in [
        "j_oh", "j_bs", "j_nb", "j_ec", "j_ic", \
        "pe_nb", "pe_ec", "pe_ic", "pe_fus", "pe_ionization", "p_rad", \
        "pi_nb", "pi_ec", "pi_ic", "pi_fus", "pi_ionization", "pi_cx", "p_ohm", "p_ei", \
        "torque_nb", "torque_in", "se_nb", "se_ionization", "si_nb", "si_ionization", \
        "q", "psipol", \
        "density_beam", "wbeam", "density_alpha", "walpha", \
        "chie",
        "chii" ]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)


    density_beam = instate["density_beam"]
    density_alpha = instate["density_alpha"]
    wbeam = instate["wbeam"]

    rho = ps["rho"]

    #------------------------------------------------------------------
    # density

    a=0; b=0; c=0; d=0
    for k in range(n_imp):
        b = b+f_imp[k]*z_imp[k]
        d = d+f_imp[k]*z_imp[k]*z_imp[k]
    for k in range(n_ion):
        a = a+f_ion[k]*z_ion[k]
        c = c+f_ion[k]*z_ion[k]*z_ion[k]

    density_ion = {}
    for k in range(n_ion):
        density_ion[k] = zeros(nrho)

    density_imp = {}
    for k in range(n_imp):
        density_imp[k] = zeros(nrho)

    density_th = zeros(nrho)

    for i in range(nrho):

        zne_adj = ne[i]
        zzne_adj = ne[i]*zeff[i]

        # depletion due to beam ions

        zne_adj = zne_adj - 1.0*density_beam[i]
        zzne_adj = zzne_adj - 1.0**2*density_beam[i]

        # effective main ion and impurity densities

        nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
        nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

        for k in range(n_ion):
            density_ion[k][i] = f_ion[k]*nion
        for k in range(n_imp):
            density_imp[k][i] = f_imp[k]*nimp

        for k in range(n_ion):
            density_th[i] = density_th[i] + density_ion[k][i]
        for k in range(n_imp):
            density_th[i] = density_th[i] + density_imp[k][i]

    ps["ns"][0] = 1.0e19*ps.node2cell(ne)

    for k in range(n_ion):
        ps["ns"][k+1] = 1.0e19*ps.node2cell(density_ion[k])

    for k in range(n_imp):
        ps["ns"][k+n_ion+1] = 1.0e19*ps.node2cell(density_imp[k])

    ps["ni"] = 1.0e19*ps.node2cell(density_th)

    #------------------------------------------------------------------
    # beam

    ps["nbeami"][0] = 1.0e19*ps.node2cell(density_beam)
    ps["eperp_beami"][0] = 2.0*20.0*ones(nrho-1)
    ps["epll_beami"][0] = 20.0*ones(nrho-1)

    tbeam = wbeam/1.602e-3/(density_beam+1.0e-6)
    ps["eperp_beami"][0] = ps.node2cell(2.0*tbeam/3.0)
    ps["epll_beami"][0] = ps.node2cell(tbeam/3.0)

    #------------------------------------------------------------------
    # temperature

    ps["Ts"][0] =ps.node2cell(te)

    for k in range(n_ion):
        ps["Ts"][k+1] = ps.node2cell(ti)
    for k in range(n_imp):
        ps["Ts"][k+n_ion+1] = ps.node2cell(ti)

    ps["Ti"] = ps.node2cell(ti)

    #------------------------------------------------------------------
    # zeff

    ps["Zeff"] = ps.node2cell(zeff)
    ps["Zeff_th"] = ps.node2cell(zeff)

    #------------------------------------------------------------------
    # rotation

    ps["omegat"] = ps.node2cell(omega)

    #--------------------------------------------------------------
    # current

    j_tot = 1.e6*instate["j_tot"]
    j_tot[0] = j_tot[1]
    ps.load_j_parallel(rho,j_tot,"rho_eq","curt",r0,b0,tot=True)

    for key in ["j_nb","j_ec","j_ic","j_bs","j_oh"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)

    j_nb  = 1.e6*instate["j_nb"]
    j_ec  = 1.e6*instate["j_ec"]
    j_ic  = 1.e6*instate["j_ic"]
    j_bs  = 1.e6*instate["j_bs"]
    j_oh  = 1.e6*instate["j_oh"]
    
    ps.load_j_parallel(rho,j_nb,"rho_nbi","curbeam",r0,b0)
    ps.load_j_parallel(rho,j_ec,"rho_ecrf","curech",r0,b0)
    ps.load_j_parallel(rho,j_ic,"rho_icrf","curich",r0,b0)
    ps.load_j_parallel(rho,j_bs,"rho","curr_bootstrap",r0,b0)
    ps.load_j_parallel(rho,j_oh,"rho","curr_ohmic",r0,b0)

    #-------------------------------------------------------------
    # heating
    
    for key in ["pe_nb","pi_nb","pe_ec","pe_ic","pi_ic","pe_fus","pi_fus"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)
    
    pe_nb  = 1.e6*instate["pe_nb" ]
    pi_nb  = 1.e6*instate["pi_nb" ]
    pe_ec  = 1.e6*instate["pe_ec" ]
    pe_ic  = 1.e6*instate["pe_ic" ]
    pi_ic  = 1.e6*instate["pi_ic" ]
    pe_fus = 1.e6*instate["pe_fus"]
    pi_fus = 1.e6*instate["pi_fus"]
    
    ps.load_vol_profile (rho,pe_nb,"rho_nbi","pbe")
    ps.load_vol_profile (rho,pi_nb,"rho_nbi","pbi")
    ps.load_vol_profile (rho,pe_ec,"rho_ecrf","peech")
    ps.load_vol_profile (rho,pe_ic,"rho_icrf","picrf_totals",k=0)
    ps.load_vol_profile (rho,pi_ic,"rho_icrf","picth")

    #------------------------------------------------------------------
    # temp

    # ps["sc0"][:] = 2.0e21
    # ps["n0norm"][:] = 1.0e-10 
    # ps["T0sc0"][:] = 0.01 
    # ps["sc0_to_sgas"][:] = 1
