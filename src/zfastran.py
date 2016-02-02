#!/usr/bin/env python

"""
 ----------------------------------------------------------------------
 utils for fastran IO
 ----------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *
import netCDF4

#--- zcode libraries
from Namelist import Namelist
from zinterp import zinterp
import zfdat, zefitutil
import zplasmastate

#---------------------------------------------------------------------- 
# variable mapping
#

"""
_fastranData = {
    "ip"        :["ip"      ,1.0] ,
    "bcentr"    :["b0"      ,1.0] ,
    "rmajor"    :["rmajor"  ,1.0] ,
    "aminor"    :["aminor"  ,1.0] ,
    "elongb"    :["kappa"   ,1.0] ,
    "trianb"    :["delta"   ,1.0] ,
    "rho"       :["rho"     ,1.0] ,
    "te"        :["te"      ,1.0] ,
    "ti"        :["ti"      ,1.0] ,
    "ene"       :["ne"      ,1.0] ,
    "amain"     :["amain"   ,1.0] ,
    "zmain"     :["zmain"   ,1.0] ,
    "aimp"      :["aimp"    ,1.0] ,
    "zimp"      :["zimp"    ,1.0] ,
    "zeff"      :["zeff"    ,1.0] ,
    "omega"     :["omega"   ,1.0] ,
    "q"         :["q"       ,1.0] ,
    "fp"        :["psipol"  ,1.0] ,
    "curpar"    :["j_tot"   ,1.0] ,
    "curbeam"   :["j_nb"    ,1.0] ,
    "curboot"   :["j_bs"    ,1.0] ,
    "currf"     :["j_rf"    ,1.0] ,
    "enbeam"    :["density_beam",1.0] ,
    "wbeam"     :["wbeam"   ,1.0] ,
    "qbeame"    :["pe_nb"   ,1.0] ,
    "qbeami"    :["pi_nb"   ,1.0] ,
    "enalp"     :["density_alpha"  ,1.0] ,
    "walp"      :["walpha"  ,1.0] ,
    "qtfuse"    :["pe_fus"  ,1.0] ,
    "qtfusi"    :["pi_fus"  ,1.0] ,
    "qrfe"      :["pe_rf"   ,1.0] ,
    "qrfi"      :["pi_rf"   ,1.0] ,
    "qrad"      :["p_rad"   ,1.0] ,
    "qohm"      :["p_ohm"   ,1.0] ,
    "qione"     :["pe_ionization",1.0] ,
    "qioni"     :["pi_ionization",1.0] ,
    "qcx"       :["pi_cx"   ,1.0] ,
    "storqueb"  :["torque_nb",1.0] ,
    "storque"   :["torque_in",1.0] ,
    "sion"      :["si_nb"   ,1.0] ,
    "chie"      :["chie"    ,1.0] ,
    "chii"      :["chii"    ,1.0] ,
}
"""

# ====================================================================== 
# io
#

def io_write_input(f_state,f_eqdsk,f_bc,rdir='.'):

    #------------------------------------------------------------------- 
    # read geqdsk and plasma state file

    geq = zefitutil.readg(f_eqdsk) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    ps = zplasmastate.zplasmastate('ips',1)
    ps.read(f_state)

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27

    z_ion = [round(x) for x in ps["qatom_S"][1:]/ps_xe ]
    a_ion = [round(x) for x in ps["m_S"][1:]/ps_mp ]
    n_ion = len(z_ion)
    n_imp = len(ps["m_SIMPI"])

    nrho  = len(ps["rho"])
    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]*1.0e-19
    nith  = ps["ni"][:]*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    zeff  = ps["Zeff"][:]
    omega = ps["omegat"][:]
    ne    = ps.cell2node_bdry(ne)
    nith  = ps.cell2node_bdry(nith)
    te    = ps.cell2node_bdry(te)
    ti    = ps.cell2node_bdry(ti)
    zeff  = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)

    q  = zeros(nrho)
    fp = zeros(nrho)

    j_tot = 1.e-6*ps.dump_j_parallel(rho,"rho_eq","curt",r0,b0,tot=True)

    j_nb = 1.e-6*ps.dump_j_parallel(rho,"rho_nbi","curbeam",r0,b0)
    j_rf = 1.e-6*(ps.dump_j_parallel(rho,"rho_ecrf","curech",r0,b0)+ps.dump_j_parallel(rho,"rho_icrf","curich",r0,b0))
    j_bs  = zeros(nrho)

    density_beam = ps.dump_profile(rho,"rho_nbi","nbeami",k=0)*1.e-19
    wbeam = ps.dump_profile(rho,"rho_nbi","eperp_beami",k=0) \
          + ps.dump_profile(rho,"rho_nbi","epll_beami" ,k=0)
    wbeam = density_beam*wbeam*1.602e-3 #MJ/m**3

    pe_nb  = ps.dump_vol_profile(rho,"rho_nbi","pbe" )*1.e-6
    pi_nb  = ps.dump_vol_profile(rho,"rho_nbi","pbi" )*1.e-6
    pth_nb = ps.dump_vol_profile(rho,"rho_nbi","pbth")*1.e-6

    density_alpha = ps.dump_profile(rho,"rho_fus","nfusi",k=0)*1.e-19
    walpha = ps.dump_profile(rho,"rho_fus","eperp_fusi",k=0) \
         + ps.dump_profile(rho,"rho_fus","epll_fusi",k=0)
    walpha = density_alpha*walpha*1.602e-3 #MJ/m**3

    pe_fus  = ps.dump_vol_profile(rho,"rho_fus","pfuse" )*1.e-6
    pi_fus  = ps.dump_vol_profile(rho,"rho_fus","pfusi" )*1.e-6
    pth_fus = ps.dump_vol_profile(rho,"rho_fus","pfusth")*1.e-6

    pe_rf  = ps.dump_vol_profile(rho,"rho_ecrf","peech") \
           + ps.dump_vol_profile(rho,"rho_icrf","pmine") \
           + ps.dump_vol_profile(rho,"rho_icrf","picrf_totals",k=0)
    pe_rf *= 1.e-6

    pi_rf  = ps.dump_profile(rho,"rho_icrf","pmini") \
           + ps.dump_profile(rho,"rho_icrf","picth")
    pi_rf *= 1.e-6

    p_rad= zeros(nrho)
    p_ohm= zeros(nrho)
    pe_ionization= zeros(nrho)
    pi_ionization= zeros(nrho)
    pi_cx= zeros(nrho)
    torque_nb= zeros(nrho)
    torque_in= zeros(nrho)
    si_nb= zeros(nrho)
    chie= zeros(nrho)
    chii= zeros(nrho)

    te[-1] = abs(te[-1])
    ti[-1] = abs(ti[-1])

    #------------------------------------------------------------------ 
    # post read plasma state

    ni = {}
    for k in range(n_ion):
        ns = ps["ns"][k+1,:]*1.0e-19
        ni[k] = ps.cell2node(ns)

    amain = zeros(nrho)
    zmain = zeros(nrho)
    aimp = zeros(nrho)
    zimp = zeros(nrho)

    for i in range(nrho):
       amain[i] = 0.0
       zmain[i] = 0.0
       sum = 0.0
       for k in range(n_ion-n_imp):
           amain[i] += a_ion[k]*ni[k][i]
           zmain[i] += z_ion[k]*ni[k][i]
           sum += ni[k][i]
       amain[i] /= sum
       zmain[i] /= sum

    for i in range(nrho):
       aimp[i] = 0.0
       zimp[i] = 0.0
       sum = 0.0
       for k in range(n_ion-n_imp,n_ion):
           aimp[i] += a_ion[k]*ni[k][i]
           zimp[i] += z_ion[k]*ni[k][i]
           sum += ni[k][i]
       aimp[i] /= sum
       zimp[i] /= sum

    time = 0.0

    #------------------------------------------------------------------ 
    # metrics

    psi     = ps["psipol"][:]/ps["psipol"][-1]  # equi-drho grid
    rho     = sqrt(ps["phit"][:]/ps["phit"][-1])
    rhob    = (ps["phit"][-1]/pi/b0)**0.5 
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

    #------------------------------------------------------------------ 
    # write inprof

    f = open(os.path.join(rdir,"inprof"),"w")

    f.write(" # generated by zfastran.py Ver+Oct2009\n")

    zfdat.write_f('inflag'  , [1.0],'',f)
    zfdat.write_f('time0'   , [time*1.0e-3],'s',f)
    zfdat.write_f('ip'      , [geq['cpasma']*1.0e-6],'',f )
    zfdat.write_f('bcentr'  , [geq["bcentr"]],'',f )
    zfdat.write_f('rmajor'  , [r0],'',f )
    zfdat.write_f('aminor'  , [rminor[-1]],'',f )
    zfdat.write_f('elongb'  , [kappa[-1]],'',f )
    zfdat.write_f('trianb'  , [delta[-1]],'',f )
    zfdat.write_f('rho'     , rho,'',f)
    zfdat.write_f('amain'   , amain,'',f)
    zfdat.write_f('zmain'   , zmain,'',f)
    zfdat.write_f('aimp'    , aimp ,'',f)
    zfdat.write_f('zimp'    , zimp ,'',f)
    zfdat.write_f('ene'     , ne,'',f)
    zfdat.write_f('te'      , te,'',f)
    zfdat.write_f('ti'      , ti,'',f)
    zfdat.write_f('zeff'    , zeff,'',f)
    zfdat.write_f('omega'   , omega,'',f)
    zfdat.write_f('q'       , q,'',f)
    zfdat.write_f('fp'      , fp,'',f)
    zfdat.write_f("curpar"  , j_tot,'',f)
    zfdat.write_f("curbeam" , j_nb,'',f)
    zfdat.write_f("currf"   , j_rf,'',f)
    zfdat.write_f("curboot" , j_bs,'',f)
    zfdat.write_f("enbeam"  , density_beam, '', f)
    zfdat.write_f("wbeam"   , wbeam,'', f)  
    zfdat.write_f("qbeame"  , pe_nb,'', f)  
    zfdat.write_f("qbeami"  , pi_nb,'', f)  
    zfdat.write_f("enalp"   , density_alpha, '', f)
    zfdat.write_f("walp"    , walpha, '', f)
    zfdat.write_f("qtfuse"  , pe_fus, '', f)
    zfdat.write_f("qtfusi"  , pi_fus, '', f)
    zfdat.write_f("qrfe"    , pe_rf,'', f)
    zfdat.write_f("qrfi"    , pi_rf,'', f)
    zfdat.write_f("qrad"    , p_rad,'', f)
    zfdat.write_f("qohm"    , p_ohm,'', f)
    zfdat.write_f("qione"   , pe_ionization,'', f)
    zfdat.write_f("qioni"   , pi_ionization,'', f)
    zfdat.write_f("qcx"     , pi_cx,'', f)
    zfdat.write_f("storqueb", torque_nb,'', f)
    zfdat.write_f("storque" , torque_in,'', f)
    zfdat.write_f("sion"    , si_nb,'', f)
    zfdat.write_f("chie"    , chie,'', f)
    zfdat.write_f("chii"    , chii,'', f)

    f.close()

    #------------------------------------------------------------------ 
    # write inmetric

    f = open(os.path.join(rdir,"inmetric"),"w")

    zfdat.write_f('Ip'     , [ip]        , '', f)
    zfdat.write_f('bcentr' , [b0]        , '', f)
    zfdat.write_f('rmajor' , [r0]        , '', f)
    zfdat.write_f('aminor' , [rminor[-1]], '', f)
    zfdat.write_f('kappa'  , [kappa [-1]], '', f)
    zfdat.write_f('delta'  , [delta [-1]], '', f)
    zfdat.write_i('nrho'   , [nrho]      , '', f)
    zfdat.write_f('rhob'   , [rhob]      , '', f)
    zfdat.write_f('rho'    , rho         , '', f)
    zfdat.write_f('volp'   , volp        , '', f)
    zfdat.write_f('ipol'   , ipol        , '', f)
    zfdat.write_f('g11'    , g11         , '', f)
    zfdat.write_f('g22'    , g22         , '', f)
    zfdat.write_f('g33'    , g33         , '', f)
    zfdat.write_f('gradrho', gradrho     , '', f)
    zfdat.write_f('area'   , area        , '', f)
    zfdat.write_f('a'      , rminor      , '', f)
    zfdat.write_f('rtor'   , rmajor      , '', f)
    zfdat.write_f('shift'  , shift       , '', f)
    zfdat.write_f('elong'  , kappa       , '', f)
    zfdat.write_f('triag'  , delta       , '', f)
    zfdat.write_f('pmhd'   , pmhd        , '', f)
    zfdat.write_f('qmhd'   , qmhd        , '', f)
    zfdat.write_f('er'     , zeros(nrho) , '', f) #<======
    zfdat.write_f('nc1'    , nc1         , '', f) #<======
    zfdat.write_f('hfac1'  , hfac1       , '', f) #<======
    zfdat.write_f('hfac2'  , hfac2       , '', f) #<======

    f.close()

    try:
        inbc = Namelist(f_bc)
        wtop = inbc["ped"]["wtop"][0]
        itop = where(rho >= 1.0-wtop)[0][0]
        infastran = Namelist("infastran")
        infastran["infastran"]["nrho_bdry"] = [itop-1]
        infastran["infastran"]["nrho_transport"] = [itop-1]
        infastran.write("infastran")
    except:
        pass

def io_update_state(
        f_state,f_eqdsk,f_fastran,time='',relax=1.0):

    #------------------------------------------------------------------ 
    # read geqdsk and plasma state

    ps = zplasmastate.zplasmastate('ips',1)
    ps.read(f_state)

    geq = zefitutil.readg(f_eqdsk) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    #------------------------------------------------------------------ 
    # read fastran

    fastran = netCDF4.Dataset(f_fastran,'r',format='NETCDF4')
    rho = fastran.variables["rho"][:]
    nrho = len(rho)

    if nrho != len(ps["rho"]):
        raise Exception('nrho differ, fastran: %d, ps: %s'%(nrho,ps.nrho))

    te = fastran.variables["te"][-1,:]
    ti = fastran.variables["ti"][-1,:]
    omega = fastran.variables["omega"][-1,:]

    j_tot = fastran.variables["j_tot"][-1,:]*1.e6
    j_oh  = fastran.variables["j_oh"][-1,:]*1.e6
    j_bs  = fastran.variables["j_bs"][-1,:]*1.e6

    #------------------------------------------------------------------
    # update plasma state

    ps["Ts"][0] = (1.0-relax)*ps["Ts"][0] + relax*ps.node2cell(te)

    nspec_th = len(ps["Ts"])-1
    print 'nspec_th =',nspec_th
    for k in range(nspec_th):
        ps["Ts"][k+1] = (1.0-relax)*ps["Ts"][k+1] + relax*ps.node2cell(ti)

    ps["Ti"][:] = (1.0-relax)*ps["Ti"][:] + relax*ps.node2cell(ti)

    ps["omegat"][:] = (1.0-relax)*ps["omegat"][:] + relax*ps.node2cell(omega)

    ps.load_j_parallel(rho,j_tot,"rho_eq","curt",r0,b0,tot=True)
    ps.load_j_parallel(rho,j_bs,"rho","curr_bootstrap",r0,b0)
    ps.load_j_parallel(rho,j_oh,"rho","curr_ohmic",r0,b0)

    #------------------------------------------------------------------
    # write plasma state

    ps.store(f_state)

#-----------------------------------------------------------------------
# test

if __name__ == "__main__":

    pass

