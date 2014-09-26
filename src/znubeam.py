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
import zplasmastate
import zutil,zinterp

########################################################################
# io

def io_write_input_files(instate,innubeam):

    #-------------------------------------------------------------------
    # machine 

    tokamak_id = instate["instate"]["tokamak_id"][0]

    #-------------------------------------------------------------------
    # species

    nspec_ion  = instate["instate"]["n_ion"] 
    z_ion      = instate["instate"]["z_ion"]        
    a_ion      = instate["instate"]["a_ion"]      
    f_ion      = instate["instate"]["f_ion"]      
    nspec_imp  = instate["instate"]["n_imp"] 
    z_imp      = instate["instate"]["z_imp"]     
    a_imp      = instate["instate"]["a_imp"]     
    f_imp      = instate["instate"]["f_imp"]      


    #-------------------------------------------------------------------
    # "dnubeam_prof.dat"

    prof = Namelist.Namelist()

    nrho  = len  (instate["instate"]["rho"]) 
    rho   = array(instate["instate"]["rho"])        
    te    = array(instate["instate"]["te"])        
    ti    = array(instate["instate"]["ti"])        
    ne    = array(instate["instate"]["ne"])
    zeff  = array(instate["instate"]["zeff"])        
    omega = array(instate["instate"]["omega"])        
    vloop = zeros(nrho) #zeros(nrho)
    n0w   = ones(nrho)  #zeros(nrho)
    n0v   = ones(nrho)  #zeros(nrho)
    f0w   = 1.0
    t0    = ones(nrho)*0.1
    dn0out = 5e+12

    prof["dnubeam_prof"]["t0_nubeam"] = [0.00]       # [s]
    prof["dnubeam_prof"]["t1_nubeam"] = innubeam["nubeam_run"]["dt_nubeam"] # [s]

    prof["dnubeam_prof"]["shot"] = [0]
    prof["dnubeam_prof"]["tokamak_id"] = [tokamak_id]

    prof["dnubeam_prof"]["nspec_ion"] = nspec_ion
    prof["dnubeam_prof"]["z_ion"]     = z_ion
    prof["dnubeam_prof"]["a_ion"]     = a_ion
    prof["dnubeam_prof"]["f_ion"]     = f_ion
    prof["dnubeam_prof"]["nspec_imp"] = nspec_imp
    prof["dnubeam_prof"]["z_imp"]     = z_imp
    prof["dnubeam_prof"]["a_imp"]     = a_imp
    prof["dnubeam_prof"]["f_imp"]     = f_imp

    prof["dnubeam_prof"]["nrho_in"  ] = [len(instate["instate"]["rho"])] # []
    prof["dnubeam_prof"]["rho_in"   ] = rho          # input grid
    prof["dnubeam_prof"]["te_in"    ] = te           # [kev]
    prof["dnubeam_prof"]["ti_in"    ] = ti           # [kev]
    prof["dnubeam_prof"]["ne_in"    ] = ne*1.0e19    # [m^-3]
    prof["dnubeam_prof"]["zeff_in"  ] = zeff         # []
    prof["dnubeam_prof"]["omega_in" ] = omega        # [rad/sec]
    prof["dnubeam_prof"]["vloop_in" ] = vloop        # [V]
    prof["dnubeam_prof"]["n0w_in"   ] = n0w          # [m^-3]
    prof["dnubeam_prof"]["n0v_in"   ] = n0v          # [m^-3]
    prof["dnubeam_prof"]["f0w_in"   ] = [f0w]        # [number/sec]
    prof["dnubeam_prof"]["t0_in"    ] = t0           # [kev]
    prof["dnubeam_prof"]["dn0out"   ] = [dn0out]     # [number/sec]

    prof["dnubeam_difb"]["difb_0"   ] = innubeam["nbi_model"]["difb_0"   ] # [cm^2/s]
    prof["dnubeam_difb"]["difb_a"   ] = innubeam["nbi_model"]["difb_a"   ] # [cm^2/s]
    prof["dnubeam_difb"]["difb_in"  ] = innubeam["nbi_model"]["difb_in"  ] # []
    prof["dnubeam_difb"]["difb_out" ] = innubeam["nbi_model"]["difb_out" ] # []
    prof["dnubeam_difb"]["nkdifb"   ] = innubeam["nbi_model"]["nkdifb"   ] # [3]

    prof.write("dnubeam_prof.dat")

    #------------------------------------------------------------------
    # "dnubeam_nbi.dat"

    nbi = Namelist.Namelist()

    nbi["dnubeam_nbi"] = innubeam["nbi_config"]

    nbi.write("dnubeam_nbi.dat")

    #------------------------------------------------------------------
    # "dnubeam_nml.dat"

    nml = Namelist.Namelist()

    nml["nbi_init"] =innubeam["nbi_init"] # innubeam["nbi_nml"]
    nml["nbi_update"]["nltest_output"] = [1]

    nml.write("dnubeam_nml.dat")

def io_update_instate(nstep=-1,navg=5,
        f_instate='instate',
        f_outnubeam='dnubeam_prof.dat'):

    #------------------------------------------------------------------ 
    # read instate
    #

    instate = Namelist.Namelist(f_instate)
    rho_ps = instate["instate"]["rho"]
    nrho_ps = len(rho_ps)

    #------------------------------------------------------------------ 
    # outnubeam
    #
    
    if nstep < 0:

        outnubeam = Namelist.Namelist(f_outnubeam)["dnubeam_out"]
    
        nrho= outnubeam["nrho_out"][0]
        rho = array(outnubeam["rho_out"][0:nrho])
        curbdotb = array(outnubeam["curbdotb_out"][0:nrho]) #A/m^2/T
        tqbsum = array(outnubeam["tqbsum_out"][0:nrho]) # NT-M/m^3
        pbe = array(outnubeam["pbe_out"][0:nrho])  #W/m^3
        pbi = array(outnubeam["pbi_out"][0:nrho])  #W/m^3
        udenspl = array(outnubeam["udenspl_out"][0:nrho]) #J/m^3 
        udenspp = array(outnubeam["udenspp_out"][0:nrho]) #J/m^3
        bdenss = array(outnubeam["bdenss_out"][0:nrho]) #/m^3

    else:

        curbdotb_a = []
        tqbsum_a = []
        pbe_a = []
        pbi_a = []
        bdenss_a = []
        udenspl_a = []
        udenspp_a = []
        
        for k in range(nstep-navg,nstep):
        
            file = "dnubeam_out_%d.dat"%k
            dat = Namelist.Namelist(file)
            dat = dat["dnubeam_out"]
            nrho     = dat["nrho_out"][0]
            rho      = dat["rho_out"][0:nrho]
            curbdotb = array(dat["curbdotb_out"][0:nrho]) #A/m^2/T
            tqbsum   = array(dat["tqbsum_out"][0:nrho]) # NT-M/m^3
            pbe      = array(dat["pbe_out"][0:nrho])  #W/m^3
            pbi      = array(dat["pbi_out"][0:nrho])  #W/m^3
            bdenss   = array(dat["bdenss_out"][0:nrho]) #/m^3
            udenspl  = array(dat["udenspl_out"][0:nrho]) #J/m^3 
            udenspp  = array(dat["udenspp_out"][0:nrho]) #J/m^3

            curbdotb_a.append(curbdotb)
            tqbsum_a.append(tqbsum)
            pbe_a.append(pbe)
            pbi_a.append(pbi)
            bdenss_a.append(bdenss)
            udenspl_a.append(udenspl)
            udenspp_a.append(udenspp)

        curbdotb = average(curbdotb_a, axis=0) 
        tqbsum   = average(tqbsum_a, axis=0)
        pbe      = average(pbe_a, axis=0)
        pbi      = average(pbi_a, axis=0)
        bdenss   = average(bdenss_a, axis=0)
        udenspl  = average(udenspl_a, axis=0)
        udenspp  = average(udenspp_a, axis=0)

    outprof = {}
    outprof["j_nb"] = curbdotb/abs(instate["instate"]["b0"][0])*1.0e-6
    outprof["torque_nb"] = tqbsum
    outprof["pe_nb"] = pbe*1.0e-6 
    outprof["pi_nb"] = pbi*1.0e-6 
    outprof["density_beam"] = bdenss*1.0e-19 
    outprof["wbeam"] = (udenspl+udenspp)*1.0e-6

    #------------------------------------------------------------------ 
    # update instate
    #

    for p in outprof.keys():

        vec = zinterp.zinterp(rho,outprof[p],s=0)[rho_ps]
        instate["instate"][p] = vec 
        
    instate.write(f_instate)

    pass

def read_nubeam_output():
    pass

###############################################################################
# driver

def driver():
    pass

###############################################################################
# test

if __name__ == "__main__":

    instate = Namelist.Namelist("instate")
    innubeam = Namelist.Namelist("innubeam")

    write_nubeam_input(instate,innubeam)
