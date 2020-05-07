"""
 -----------------------------------------------------------------------
 utils for genray IO
 -----------------------------------------------------------------------
"""

from numpy import *
import netCDF4

from Namelist import Namelist
from efit_eqdsk import readg
from zinterp import zinterp
from plasmastate import plasmastate

def write_inputfiles(f_state, f_eqdsk, f_ingenray, MKS=True):
    #-- read plasma state file
    ps = plasmastate('ips', 1)
    ps.read(f_state)

    ps_xe = 1.6022e-19
    ps_mp = 1.6726e-27

    z_ion = ps["qatom_S"][1:]/ps_xe
    a_ion = ps["m_S"][1:]/ps_mp
    n_ion = len(z_ion)
    n_imp = len(ps["m_SIMPI"])

    nrho = len(ps["rho"])
    rho = ps["rho"][:]
    if MKS:
       ne = ps["ns"][0,:]
    else:
       ne = ps["ns"][0,:]*1.0e-19
    te = ps["Ts"][0,:]
    ti = ps["Ti"][:]
    zeff = ps["Zeff"][:]
    ne = ps.cell2node(ne)
    te = ps.cell2node(te)
    ti = ps.cell2node(ti)
    zeff = ps.cell2node(zeff)

    ni = {}
    for k in range(n_ion):
        if MKS:
            ns = ps["ns"][k+1,:]
        else:
            ns = ps["ns"][k+1,:]*1.0e-19
        ni[k] = ps.cell2node(ns)

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

    #-- genray.dat
    ingenray = Namelist(f_ingenray)
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
    if MKS:
        ingenray.write("genray.in")
    else:
        ingenray.write("genray.dat")

def update_state(f_state, f_eqdsk, imode='EC', jmulti=1.0, add=0, rho_smooth=0.0):
    #-- read genray output
    ncfile = "genray.nc"
    ncgenray = netCDF4.Dataset(ncfile, 'r', format='NETCDF4')

    rho  = ncgenray.variables['rho_bin_center'][:]
    nrho = len(rho)

    j_par = ncgenray.variables['s_cur_den_onetwo'][:]*0.01 # MA/m**2
    prf_e = ncgenray.variables['powden_e'][:]*1.0e-7 # erg/(cm**3*sec) -> MW/(m**3*sec)
    prf_i = ncgenray.variables['powden_i'][:]*1.0e-7 # erg/(cm**3*sec) -> MW/(m**3*sec)
    I_rf  = ncgenray.variables['toroidal_cur_total'][:]
    j_par = jmulti*abs(j_par)

    i = where (rho < rho_smooth)[0]
    j_par[i] = 0.0
    prf_e[i] = 0.0
    prf_i[i] = 0.0

    ncgenray.close()

    #-- update plasma state
    ps = plasmastate('ips', 1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0 = geq["rzero" ]
    b0 = abs(geq["bcentr"])

    rho_ps = ps["rho"][:]
    jp_ps = abs(1.0e6*zinterp(rho, j_par)(rho_ps))
    pe_ps = 1.0e6*zinterp(rho, prf_e)(rho_ps)
    pi_ps = 1.0e6*zinterp(rho, prf_i)(rho_ps)

    if imode=='EC':
        print('imode:EC')
        ps.load_j_parallel(rho_ps, jp_ps, "rho_ecrf", "curech", r0, b0)
        ps.load_vol_profile(rho_ps, pe_ps, "rho_ecrf", "peech")
    elif imode=='IC':
        print('imode:EC')
        ps.load_j_parallel(rho_ps, jp_ps, "rho_icrf", "curich", r0, b0, add=add)
        ps.load_vol_profile(rho_ps, pe_ps, "rho_icrf", "picrf_totals", k=0, add=add)
        ps.load_vol_profile(rho_ps, pi_ps, "rho_icrf", "picrf_totals", k=1, add=add)

    ps.store(f_state)

    print("I_genray = ",I_rf)
