#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
  fixed boundary efit solver
  JM
 -----------------------------------------------------------------------
"""

import sys,os,glob,pickle,shutil
from numpy import *
from scipy.interpolate import splrep,splev,interp1d

from Namelist import Namelist
import zutil,zefitutil,zefitdata
from zinterp import *
from zplasmastate import plasma_state_file

def fixbdry_kfile(shot,time,f_inefit,p_scale=1.0,efitdir='.'): 

    mu0 = 4.e-7*pi

    #-------------------------------------------------------------------
    # read inefit
    #

    inefit = Namelist(f_inefit,"r")

    ip = inefit["inefit"]["ip"][0] # [A] 
    r0 = inefit["inefit"]["r0"][0] # [m]
    b0_sign = inefit["inefit"]["b0"][0] # [T]
    b0 = abs(b0_sign)
    rho_in = inefit["inefit"]["rho"] # []
    press_in = inefit["inefit"]["press"] #[Pa]
    jpar_in = inefit["inefit"]["jpar"] #[A/m^2]

    nrho_in = len(rho_in)
    rho_in = array(rho_in)
    press_in = array(press_in)
    jpar_in = array(jpar_in)

    #-------------------------------------------------------------------
    # scale pressure
    #

    press_in = p_scale*press_in

    #-------------------------------------------------------------------
    # read geqdsk file and calculate metric
    #

    gfile = os.path.join(efitdir,'g%06d.%05d'%(shot,time))

    print 'efitdir = ',efitdir
    print 'gfile   = ',gfile

    geq = zefitdata.efitdata(gfile,nrho_eq=129,nth_eq=101,nrho_eq_geo=129)

    psirho = geq["psipol"]/geq["psipol"][-1]  # equi-drho grid
    dpsi = abs(geq["ssibry"]-geq["ssimag"])

    rhob = (geq["phit"][-1]/pi/abs(geq["bcentr"]))**0.5 

    nrho = len(psirho)
    rho = arange(nrho)/(nrho-1.0)
    drho = 1.0/(nrho-1.0)

    npsi = len(psirho)
    psi = arange(npsi)/(npsi-1.0)
    rhopsi_spl = interp1d(psirho,rho,kind='cubic')
    rho_eval = rhopsi_spl(psi)
    rho_eval[0]  = 0.0
    rho_eval[-1] = 1.0

    ipol = geq["g_eq"][:]/(r0*b0)
    volp = 4.0*pi*pi*rho*rhob*r0/ipol/(r0*r0*geq["gr2i"][:])
    rm2  = geq["gr2i"][:]
    qmhd = geq["q_eq"][:]

    #-------------------------------------------------------------------
    # spline profiles
    #

    ipol_spl = zinterp(rho,ipol)
    ipol_der = ipol_spl(rho,der=1)

    jpar_spl = zinterp(rho_in,jpar_in)
    jpar = jpar_spl(rho)

    rm2_spl = zinterp(rho,rm2)
    rm2_psi = rm2_spl(rho_eval)

    jtor = zeros(nrho)
    for i in range(1,nrho):
        jtor[i] = 1.0/ipol[i]*(jpar[i]
            +4.0*pi**2/mu0*r0*b0/volp[i]*rho[i]/qmhd[i]*ipol_der[i])
    jtor[0] = 1.0/ipol[0]*jpar[0]

    jtor_spl = zinterp(rho,jtor)
    jtor_psi = jtor_spl(rho_eval)

    press_spl = zinterp(rho_in,press_in,s=0) 
    press_psi = press_spl(rho_eval,der=0) 
   
    tmp = zinterp(psi,press_psi,s=0) 
    pprim_psi = tmp(psi,der=1) 
    pprim_psi = pprim_psi/dpsi

    #-------------------------------------------------------------------
    # Ip scale
    #
    
    ip_temp = 0.0
    for i in range(nrho-1):
        rho_m  = 0.5*(rho[i+1]+rho[i])
        volp_m = 0.5*(volp[i+1]+volp[i])
        drho_m = rhob*(rho[i+1]-rho[i])
        ip_temp += volp_m*jtor_spl(rho_m)*drho_m
    ip_temp /= 2.0*pi*r0
    print 'Ip = %5.3f %5.3f %5.3f'%(ip*1.0e-6,ip_temp*1.0e-6,ip/ip_temp)

    for i in range(len(jtor_psi)):
        jtor_psi[i] = ip/ip_temp*jtor_psi[i] 

    #-------------------------------------------------------------------
    # calculate p' and ff'
    #

    x0 = [0.1,0.15,0.2]
    extrapolation(psi,pprim_psi,x0,y0=None)

    ffprim_psi = zeros(npsi)
    for i in range(npsi):
       #ffprim_psi[i] = -mu0*r0/(r0**2*rm2[i])*(jtor_psi[i]+r0*pprim_psi[i]) 
        ffprim_psi[i] = -mu0*r0/(r0**2*rm2_psi[i])*(jtor_psi[i]+r0*pprim_psi[i])

    eval = extrapolation_sep(
              [psi[-4],psi[-3],psi[-2]],
              [pprim_psi[-4],pprim_psi[-3],pprim_psi[-2]])
    pprim_psi[-1] = eval

    eval = extrapolation_sep(
              [psi[-4],psi[-3],psi[-2]],
              [ffprim_psi[-4],ffprim_psi[-3],ffprim_psi[-2]])
    ffprim_psi[-1] = eval

    #-------------------------------------------------------------------
    # write kfile
    #

    kfile = Namelist()

    kfile['in1']['iconvr'] = [3]
    kfile['in1']['mxiter'] = [-1]
    kfile['in1']['nxiter'] = [51]
    kfile['in1']['error' ] = [1e-04]
    kfile['in1']['relax' ] = [0.5]
    kfile['in1']['iout'  ] = [4]
    kfile['in1']['icurrt'] = [2]
    kfile['in1']['kppfnc'] = [6]
    kfile['in1']['kppcur'] = [3]
    kfile['in1']['kfffnc'] = [6]
    kfile['in1']['kffcur'] = [3]
    kfile['in1']['ishot' ] = [shot] 
    kfile['in1']['itime' ] = [time]
    kfile['in1']['itimeu'] = [0]
    kfile['in1']['itek'  ] = [0]
    kfile['in1']['itrace'] = [1]
    kfile['in1']['plasma'] = [ip]
    kfile['in1']['btor'  ] = [b0_sign]
    kfile['in1']['rcentr'] = [r0]
    kfile['in1']['rzero' ] = [r0]
    kfile['in1']['limitr'] = inefit["inefit"]["nlim"] 
    kfile['in1']['xlim'  ] = inefit["inefit"]["rlim"]
    kfile['in1']['ylim'  ] = inefit["inefit"]["zlim"]
    kfile['in1']['nbdry' ] = inefit["inefit"]["nbdry"] 
    kfile['in1']['rbdry' ] = inefit["inefit"]["rbdry"]
    kfile['in1']['zbdry' ] = inefit["inefit"]["zbdry"]
    kfile['in1']['ierchk'] = [0]
    kfile['in1']['nextra'] = [0]
    kfile['in1']['scrape'] = [0.01]
    kfile['in1']['psibry'] = [0.0]
    kfile['in1']['cfcoil'] = [-1.0]
    kfile['in1']['enp'   ] = [1.0]
    kfile['in1']['emp'   ] = [1.0]
    kfile['in1']['fcurbd'] = [0]
    kfile['in1']['pcurbd'] = [0]
    kfile['in1']['fczero'] = 16*[1.0]
    kfile['in1']['fcsum' ] = 16*[1.0]
    kfile['in1']['ifref' ] = [2]
    kfile['in1']['prbdry'] = [press_in[-1]]
    kfile['profile_ext']['npsi_ext'] = [npsi]
    kfile['profile_ext']['psin_ext'] = psi
    kfile['profile_ext']['pprime_ext'] = -pprim_psi
    kfile['profile_ext']['ffprim_ext'] = -ffprim_psi

    kfile.write("k%06d.%05d"%(shot,time))

def fixbdry_kfile_init(shot,time,f_inefit): 

    #-------------------------------------------------------------------
    # read input 
    #

    inefit = Namelist(f_inefit,"r")
    rho_in = inefit["inefit"]["rho"]
    press_in = inefit["inefit"]["press"]
    jtor_in = inefit["inefit"]["jpar"]

    ip = inefit["inefit"]['ip'][0]  
    b0 = inefit["inefit"]['b0'][0]  
    r0  = inefit["inefit"]['r0'][0]  

    npsi = 129
    psi  = arange(npsi)/(npsi-1.0)

    #-------------------------------------------------------------------
    # write kfile
    #

    kfile = Namelist()

    kfile['in1']['iconvr'] = [3]
    kfile['in1']['mxiter'] = [-1]
    kfile['in1']['nxiter'] = [51]
    kfile['in1']['error' ] = [1e-04]
    kfile['in1']['relax' ] = [0.5]
    kfile['in1']['iout'  ] = [4]
    kfile['in1']['icurrt'] = [2]
    kfile['in1']['kppfnc'] = [6]
    kfile['in1']['kppcur'] = [3]
    kfile['in1']['kfffnc'] = [6]
    kfile['in1']['kffcur'] = [3]
    kfile['in1']['ishot' ] = [shot] 
    kfile['in1']['itime' ] = [time]
    kfile['in1']['itimeu'] = [0]
    kfile['in1']['itek'  ] = [0]
    kfile['in1']['itrace'] = [1]
    kfile['in1']['plasma'] = [ip]
    kfile['in1']['btor'  ] = [b0]
    kfile['in1']['rcentr'] = [r0]
    kfile['in1']['rzero' ] = [r0]
    kfile['in1']['limitr'] = inefit["inefit"]["nlim"] 
    kfile['in1']['xlim'  ] = inefit["inefit"]["rlim"]
    kfile['in1']['ylim'  ] = inefit["inefit"]["zlim"]
    kfile['in1']['nbdry' ] = inefit["inefit"]["nbdry"] 
    kfile['in1']['rbdry' ] = inefit["inefit"]["rbdry"]
    kfile['in1']['zbdry' ] = inefit["inefit"]["zbdry"]
    kfile['in1']['ierchk'] = [0]
    kfile['in1']['nextra'] = [0]
    kfile['in1']['scrape'] = [0.01]
    kfile['in1']['psibry'] = [0.0]
    kfile['in1']['cfcoil'] = [-1.0]
    kfile['in1']['enp'   ] = [1.0]
    kfile['in1']['emp'   ] = [1.0]
    kfile['in1']['fcurbd'] = [0]
    kfile['in1']['pcurbd'] = [0]
    kfile['in1']['fczero'] = 16*[1.0]
    kfile['in1']['fcsum' ] = 16*[1.0]
    kfile['in1']['ifref' ] = [2]
    kfile['in1']['prbdry'] = [0.0]
    kfile['profile_ext']['npsi_ext'] = [npsi]
    kfile['profile_ext']['psin_ext'] = psi
    kfile['profile_ext']['pprime_ext'] = npsi*[1.0]
    kfile['profile_ext']['ffprim_ext'] = npsi*[0.0]

    kfile.write("k%06d.%05d"%(shot,time))

########################################################################
#   UTILS

def inverse3x3(x,y):

    a = array(
      [ [x[0]**2,x[0],1.0],
        [x[1]**2,x[1],1.0],
        [x[2]**2,x[2],1.0] ] )
    a_inv = linalg.inv(a)
    return dot(a_inv,y)

def cubic(p,x0):
    return p[0]*x0**2+p[1]*x0+p[2]

def extrapolation(rho,f,x0,y0):

    spl = zinterp(rho,f,s=0) #splrep(rho,f,s=0)

    if y0 == None:
        #p = inverse3x3(x0,[spl(x) for x in x0])
        tmp = spl(x0,der=0) # splev(x0,spl,der=0)
        p = inverse3x3(x0,tmp)
        for k in range(len(rho)):
            if rho[k] < x0[0]:
                f[k] = cubic(p,rho[k])
    else:
        #p = inverse3x3(x0,[y0,spl(x0[1]),spl(x0[2])])
        tmp = spl([x0[1],x0[2]],der=0) #  splev([x0[1],x0[2]],spl,der=0)
        p = inverse3x3(x0,[y0,tmp[0],tmp[1]])
        for k in range(len(rho)):
            if rho[k] < x0[1]:
                f[k] = cubic(p,rho[k])

def extrapolation_sep(x0,y0):

    p = inverse3x3(x0,y0)
    return cubic(p,1.0)     

########################################################################
#  STATE IO

def io_input_from_instate(f_instate):

    instate = Namelist(f_instate)
    instate = instate['instate']

    nrho  = instate['nrho'][0]
    rho   = array(instate['rho']) 

    ne    = array(instate['ne'])
    te    = array(instate['te'])
    ti    = array(instate['ti'])
    wbeam = array(instate['wbeam'])
    walp  = array(instate['walpha'])
    zeff  = array(instate['zeff'])

    z_ion = 0.0
    a_ion = 0.0
    for k in range(instate["n_ion"][0]):
        z_ion += instate['z_ion'][k]*instate['f_ion'][k]
        a_ion += instate['a_ion'][k]*instate['f_ion'][k]

    z_imp = 0.0
    a_imp = 0.0
    for k in range(instate["n_imp"][0]):
        z_imp += instate['z_imp'][k]*instate['f_imp'][k]
        a_imp += instate['a_imp'][k]*instate['f_imp'][k]

    zmain = z_ion*ones(nrho)
    zimp  = z_imp*ones(nrho) 

    nalpha = array(instate['density_alpha'])
    nbfast = array(instate['density_beam'])
    nhe = zeros(nrho)

    jpar  = array(instate['j_tot'])*1.0e6

    tmp = zmain*zimp*(zimp-zmain)
    ni = (zimp**2*(ne-nbfast-2.0*(nalpha+nhe))
         -zimp*(zeff*ne-nbfast-4.0*(nalpha+nhe))
             )/tmp;
    nz = (ne*(zeff-1.0)-2.0*(nalpha+nhe))/tmp;

    pmhd = 1.602e3*(ne*te +(ni+nz)*ti) \
          +2.0/3.0*1.0e6*(wbeam+walp)

    inefit = Namelist()
    
    inefit["inefit"]["ip"   ] = [instate["ip"][0]*1.0e6]
    inefit["inefit"]["r0"   ] = instate["r0"]
    inefit["inefit"]["b0"   ] = instate["b0"]
    inefit["inefit"]["nrho" ] = [nrho]
    inefit["inefit"]["rho"  ] = rho
    inefit["inefit"]["press"] = pmhd
    inefit["inefit"]["jpar" ] = jpar
    inefit["inefit"]["nlim" ] = instate["nlim" ]
    inefit["inefit"]["rlim" ] = instate["rlim" ]
    inefit["inefit"]["zlim" ] = instate["zlim" ]
    inefit["inefit"]["nbdry"] = instate["nbdry"]
    inefit["inefit"]["rbdry"] = instate["rbdry"]
    inefit["inefit"]["zbdry"] = instate["zbdry"]

    inefit.write("inefit")

def io_input_from_state(f_geqdsk,f_ps,f_inbc):

    # read inbc

    inbc = Namelist(f_inbc)["inbc"]

    # read geqdsk

    geq = zefitutil.readg(f_geqdsk) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    # read plasma state

    ps = plasma_state_file(f_ps,r0=r0,b0=b0,ip=ip)

    nrho  = ps.nrho

    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]*1.0e-19
    ni    = ps["ni"][:]*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    ne    = ps.cell2node(ne)
    ni    = ps.cell2node(ni)
    te    = ps.cell2node(te)
    ti    = ps.cell2node(ti)

    density_beam = ps.dump_profile(rho,"nbeami",k=0)*1.e-19
    wbeam = ps.dump_profile(rho,"eperp_beami",k=0) \
          + ps.dump_profile(rho,"epll_beami" ,k=0)
    wbeam = density_beam*wbeam

    density_alpha = ps.dump_profile(rho,"nfusi",k=0)*1.e-19
    walpha = ps.dump_profile(rho,"eperp_fusi",k=0) \
           + ps.dump_profile(rho,"epll_fusi",k=0)
    walpha = density_alpha*walpha

    pmhd = 1.602e3*(ne*te+ni*ti)+2.0/3.0*1.602e3*(wbeam+walpha)

    #print 1.602e3*(ne*te+ni*ti)
    #print 2.0/3.0*(wbeam+walpha)

    jpar = ps.dump_j_parallel()

    rlim = ps["rlim"][:]
    zlim = ps["zlim"][:]
    nlim = len(rlim)

    ps.close()

    # namelist for efit

    inefit = Namelist()
    inefit["inefit"]["ip"   ] = [inbc["ip"][0]*1.0e6]
    inefit["inefit"]["r0"   ] = inbc["r0"]
    inefit["inefit"]["b0"   ] = inbc["b0"]
    inefit["inefit"]["nrho" ] = [nrho]
    inefit["inefit"]["rho"  ] = rho
    inefit["inefit"]["press"] = pmhd
    inefit["inefit"]["jpar" ] = jpar
    inefit["inefit"]["nlim" ] = [nlim]
    inefit["inefit"]["rlim" ] = rlim
    inefit["inefit"]["zlim" ] = zlim
    inefit["inefit"]["nbdry"] = inbc["nbdry"]
    inefit["inefit"]["rbdry"] = inbc["rbdry"]
    inefit["inefit"]["zbdry"] = inbc["zbdry"]
    inefit.write("inefit")


########################################################################
#  STAND ALONE   

if __name__ == "__main__":

    #---------------------------------------------------
    # run time options

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("--shot",
        action="store",type="int",dest="shot")
    parser.add_option("--time",
        action="store",type="int",dest="time")
    parser.add_option("--efitdir",
        action="store",type="string",dest="efitdir",default=".")
    parser.add_option("--niter",
        action="store",type="int",dest="niter",default=5)
    (options,args) = parser.parse_args(sys.argv[1:])

    shot = options.shot 
    time = options.time
    efitdir = options.efitdir
    niter = options.niter

    runefit = 'efitd90 129 129'

    f_inefit='inefit'
            
    for k in range(niter):

        zutil.outscreen("efit iteration = %d/%d"%(k,niter))

        print "generate kfile"
        if k == 0:
            fixbdry_kfile_init(shot,time,f_inefit)
        else:
            fixbdry_kfile(shot,time,f_inefit)
    
        print "run efit"
        kfilename = "k%06d.%05d"%(shot,time)
        f = open("efparm","w")
        f.write("2\n1\n%s\n"%kfilename)
        f.close()
        os.system(runefit+" <efparm"+" >& elog%05d_%d"%(time,k)) 
