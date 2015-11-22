#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 utils for genray IO
 JM
 -----------------------------------------------------------------------
"""

import os,sys,shutil
from numpy import *
import netCDF4

#-------------------
#--- zcode libraries
from Namelist import Namelist
import zefitutil
from zinterp import zinterp
from zplasmastate import plasma_state_file

def io_write_inputfiles(f_state,f_eqdsk,f_ingenray):

    # =================================================================
    # read plasma state file
    
    ps  = plasma_state_file(f_state)

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27

    z_ion = ps["qatom_S"][1:]/ps_xe
    a_ion = ps["m_S"][1:]/ps_mp
    n_ion = len(z_ion)
    n_imp = len(ps.data.dimensions["dim_nspec_imp0"])

    nrho  = ps.nrho
    rho   = ps["rho"][:]
### ne    = ps["ns"][0,:]*1.0e-19
    ne    = ps["ns"][0,:]
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    zeff  = ps["Zeff"][:]
    ne    = ps.cell2node(ne)
    te    = ps.cell2node(te)
    ti    = ps.cell2node(ti)
    zeff  = ps.cell2node(zeff)

    ni = {}
    for k in range(n_ion):
        ns = ps["ns"][k+1,:]*1.0e-19
        ni[k] = ps.cell2node(ns)

    # =================================================================
    # genray.dat

    ingenray = Namelist(f_ingenray)

    nbulk = 1 + n_ion
    prof_den = []
    prof_tmp = []
    prof_den.append(ne)
    prof_tmp.append(te)
    for k in range(n_ion):
        prof_den.append(ni[k])
        prof_tmp.append(ti)
    prof_den = ravel(array(prof_den).transpose())
    prof_tmp = ravel(array(prof_tmp).transpose())

    ingenray["TOKAMAK"]["EQDSKIN"] = [f_eqdsk]

    ingenray["PLASMA"]["NDENS"] = [nrho]
    ingenray["PLASMA"]["NBULK"] = [nbulk] 
    ingenray["PLASMA"]["IZEFF"] = [2] 
    ingenray["PLASMA"]["IDENS"] = [1] 
    ingenray["PLASMA"]["TEMP_SCALE"] = nbulk*[1.0] 
    ingenray["PLASMA"]["DEN_SCALE"] = nbulk*[1.0] 

    ingenray["DENTAB"]["PROF"] = prof_den
    ingenray["TEMTAB"]["PROF"] = prof_tmp
    ingenray["ZEFTAB"]["ZEFF1"] = zeff

#   ingenray.write("genray.dat")
    ingenray.write("genray.in")

def io_update_state(f_state,f_eqdsk,imode='IC'):

    # read genray output

    ncfile = "genray.nc" 
    ncgenray = netCDF4.Dataset(ncfile,'r',format='NETCDF4')

    rho   = ncgenray.variables['rho_bin_center'][:]

    j_par = ncgenray.variables['s_cur_den_onetwo'][:]*0.01 # MA/m**2
    prf_e = ncgenray.variables['powden_e'][:]*1.0e-7  # erg/(cm**3*sec) -> MW/(m**3*sec)
    prf_i = ncgenray.variables['powden_i'][:]*1.0e-7  # erg/(cm**3*sec) -> MW/(m**3*sec)
    I_rf  = ncgenray.variables['toroidal_cur_total'][:]
    j_par = abs(j_par)
   #j_par = -j_par

    ncgenray.close()

    # update plasma state

    geq = zefitutil.readg(f_eqdsk) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']
    ps  = plasma_state_file(f_state,r0=r0,b0=b0,ip=ip)

    rho_ps = ps["rho"][:]
    jp_ps = 1.0e6*zinterp(rho,j_par)(rho_ps)
    pe_ps = 1.0e6*zinterp(rho,prf_e)(rho_ps)

#   ps.load_j_parallel_CD(rho_ps,jp_ps,"ec")
#   ps.load_profile(rho_ps,pe_ps,"peech","vol")

    ps.load_j_parallel_CD(rho_ps,jp_ps,"ic")
    ps.load_profile(rho_ps,pe_ps,"pmine","vol")
#   ps.load_profile(rho_ps,pi_ps,"pmini","vol")


#   ps.load_j_parallel_CD(rho_ps,jp_ps,"ec")
#   ps.load_profile(rho_ps,pe_ps,"peech","vol")


    ps.close()

    print "I_genray = ",I_rf
