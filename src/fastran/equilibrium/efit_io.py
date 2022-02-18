"""
 -----------------------------------------------------------------------
  fixed boundary efit solver
 -----------------------------------------------------------------------
"""
import shutil
from numpy import *
from scipy.interpolate import interp1d
import time as timer
from Namelist import Namelist
from fastran.util.zinterp import zinterp
from fastran.plasmastate.plasmastate import plasmastate

########################################################################
#   compose kfile

def fixbdry_kfile(shot, time, inefit, inmetric, relax=0, topology='', error=1.0e-4):
    mu0 = 4.e-7*pi

    #-------------------------------------------------------------------
    # from inefit
    ip = inefit["inefit"]["ip"][0] # [A]
    r0 = inefit["inefit"]["r0"][0] # [m]
    b0_sign = inefit["inefit"]["b0"][0] # [T]
    b0 = abs(b0_sign)

    rho_in = array(inefit["inefit"]["rho"]) # []
    press_in = array(inefit["inefit"]["press"]) #[Pa]
    jpar_in = array(inefit["inefit"]["jpar"]) #[A/m^2]

    #-------------------------------------------------------------------
    # from inmetric
    for key in inmetric.keys():
        inmetric[key] = array(inmetric[key])

    psirho = inmetric["psi"]/inmetric["psi"][-1]
    dpsi = abs(inmetric["psi"][0]-inmetric["psi"][-1])

    rhob = inmetric["rhob"][0]
    nrho = inmetric["nrho"][0]
    rho = inmetric["rho"]
    drho = 1.0/(nrho-1.0)

    npsi = nrho
    psi = arange(npsi)/(npsi-1.0)

    rhopsi_spl = interp1d(psirho, rho, kind='cubic')
    rho_eval = rhopsi_spl(psi)
    rho_eval[0]  = 0.0
    rho_eval[-1] = 1.0

    ipol = inmetric["ipol"]
    vol  = inmetric["vol"]
    rm2  = inmetric["gr2i"]

    #-------------------------------------------------------------------
    # spline profiles
    ipol_spl = zinterp(rho, ipol)
    ipol_der = ipol_spl(rho, der=1)

    jpar_spl = zinterp(rho_in, jpar_in)
    jpar = jpar_spl(rho)

    rm2_spl = zinterp(rho, rm2)
    rm2_psi = rm2_spl(rho_eval)

    curt = zeros(nrho)
    curt[0] = 0.0
    for i in range(1,nrho):
        dV = (vol[i]-vol[i-1])
        ipol_m = 0.5*(ipol[i]+ipol[i-1])
        jpar_m = 0.5*(jpar[i]+jpar[i-1])
        curt[i] = (curt[i-1]/ipol[i-1] + jpar_m/ipol_m**2/(2.0*pi*r0) * dV)*ipol[i]

    jtor = zeros(nrho+1)
    rho_cell = zeros(nrho+1)
    for i in range(1,nrho):
        dV = (vol[i]-vol[i-1])
        jtor[i] = (curt[i]-curt[i-1])/dV*2.0*pi*r0
        rho_cell[i] = 0.5*(rho[i]+rho[i-1])
    jtor[0] = jtor[1]
    jtor[-1] = 2.0*jtor[-2]-jtor[-3]
    rho_cell[0] = 0.0
    rho_cell[-1] = 1.0

    jtor_spl = zinterp(rho_cell, jtor)
    jtor_psi = jtor_spl(rho_eval)

    press_spl = zinterp(rho_in, press_in, s=0)
    press_psi = press_spl(rho_eval, der=0)

    tmp = zinterp(psi, press_psi, s=0)
    pprim_psi = tmp(psi, der=1)
    pprim_psi = pprim_psi/dpsi

    ipol_psi = ipol_spl(rho_eval, der=0) #------

    #-------------------------------------------------------------------
    # Ip scale
    ip_temp = 0.0
    for i in range(nrho-1):
        dV = (vol[i+1]-vol[i])
        rho_m  = 0.5*(rho[i+1]+rho[i])
        ip_temp += dV*jtor_spl(rho_m)
    ip_temp /= 2.0*pi*r0
    print ('ip, ip cal, ip/(ip cal) = %5.3f %5.3f %5.3f'%(ip*1.0e-6, ip_temp*1.0e-6, ip/ip_temp))

    for i in range(len(jtor_psi)):
        jtor_psi[i] = ip/ip_temp*jtor_psi[i]

    #-------------------------------------------------------------------
    # calculate p' and ff'
    x0 = [0.1,0.15,0.2]
    extrapolation(psi, pprim_psi, x0, y0=None)

    ffprim_psi = zeros(npsi)
    for i in range(npsi):
        ffprim_psi[i] = -mu0*r0/(r0**2*rm2_psi[i])*(jtor_psi[i]+r0*pprim_psi[i])

    eval = extrapolation_sep(
              [psi[-4], psi[-3], psi[-2]],
              [pprim_psi[-4], pprim_psi[-3], pprim_psi[-2]])
    pprim_psi[-1] = eval

    eval = extrapolation_sep(
              [psi[-4], psi[-3], psi[-2]],
              [ffprim_psi[-4], ffprim_psi[-3], ffprim_psi[-2]])
    ffprim_psi[-1] = eval

    #-------------------------------------------------------------------
    # x point

    rbdry = inefit["inefit"]["rbdry"]
    zbdry = inefit["inefit"]["zbdry"]

    k_x1 = argmax(zbdry)
    k_x2 = argmin(zbdry)

    #-------------------------------------------------------------------
    # write kfile

    if relax:
       kfile_old = Namelist("k%06d.%05d"%(shot, time))
       infreegs_old = Namelist("infreegs")

    kfile = Namelist()

    kfile['in1']['iconvr'] = [3]
    kfile['in1']['mxiter'] = [-1]
    kfile['in1']['nxiter'] = [51]
    kfile['in1']['error' ] = [error]
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
    kfile['in1']['cfcoil'] = [-10.0] #[-1.0]
    kfile['in1']['enp'   ] = [1.0]
    kfile['in1']['emp'   ] = [1.0]
    kfile['in1']['fcurbd'] = [0]
    kfile['in1']['pcurbd'] = [0]
    kfile['in1']['fczero'] = 16*[1.0]
    kfile['in1']['fcsum' ] = 16*[1.0]
    kfile['in1']['ifref' ] = [2]
    kfile['in1']['prbdry'] = [press_in[-1]]
    kfile['in1']['kbetapr'] = [1]

    if topology in ['LSN', 'USN', 'DN'] :
       kfile['in1']['fwtbdry'] = inefit["inefit"]["nbdry"][0] * [1.0]
    if topology in ['USN', 'DN'] :
        print ('weigting upper xpoint:', list(range(k_x1-1, k_x1+3)))
        for k in range(k_x1-1, k_x1+3):
            if k<=0:
               _k =  inefit["inefit"]["nbdry"][0] + k
            else:
               _k = k
            kfile['in1']['fwtbdry(%d)'%_k] = [10.]
    if topology in ['LSN', 'DN']:
        print ('weigting lower xpoint:', list(range(k_x2-1, k_x2+3)))
        for k in range(k_x2-1, k_x2+3):
            if k<=0:
               _k =  inefit["inefit"]["nbdry"][0] + k
            else:
               _k = k
            kfile['in1']['fwtbdry(%d)'%_k] = [10.]

    kfile['profile_ext']['npsi_ext'] = [npsi]
    kfile['profile_ext']['psin_ext'] = psi
    if relax:
        kfile['profile_ext']['pprime_ext'] = 0.5*(-pprim_psi + kfile_old['profile_ext']['pprime_ext'])
        kfile['profile_ext']['ffprim_ext'] = 0.5*(-ffprim_psi + kfile_old['profile_ext']['ffprim_ext'])
    else:
        kfile['profile_ext']['pprime_ext'] = -pprim_psi
        kfile['profile_ext']['ffprim_ext'] = -ffprim_psi

    kfile.write("k%06d.%05d"%(shot, time))

    infreegs = Namelist()
    infreegs['profile_ext']['npsi_ext'] = [npsi]
    infreegs['profile_ext']['psin_ext'] = psi
    if relax:
        infreegs['profile_ext']['pprime_ext'] = 0.5*(-pprim_psi + infreegs_old['profile_ext']['pprime_ext'])
        infreegs['profile_ext']['ffprim_ext'] = 0.5*(-ffprim_psi + infreegs_old['profile_ext']['ffprim_ext'])
        infreegs['profile_ext']['p_ext'] = 0.5*(press_psi + infreegs_old['profile_ext']['p_ext'])
        infreegs['profile_ext']['f_ext'] = 0.5*(ipol_psi*r0*b0 + infreegs_old['profile_ext']['f_ext'])
    else:
        infreegs['profile_ext']['pprime_ext'] = -pprim_psi
        infreegs['profile_ext']['ffprim_ext'] = -ffprim_psi
        infreegs['profile_ext']['p_ext'] = press_psi
        infreegs['profile_ext']['f_ext'] = ipol_psi*r0*b0
    infreegs.write("infreegs")

    #-------------------------------------------------------------------
    # convergence check

    iconv = 0
    max_error = 0.01
    diff_p = 0.0
    diff_p_sum = 0.0
    diff_f = 0.0
    diff_f_sum = 0.0
    if relax:
       for k in range(nrho):
           diff_p += abs(kfile['profile_ext']['pprime_ext'][k] - kfile_old['profile_ext']['pprime_ext'][k] )
           diff_p_sum += abs(kfile['profile_ext']['pprime_ext'][k])
           diff_f += abs(kfile['profile_ext']['ffprim_ext'][k] - kfile_old['profile_ext']['ffprim_ext'][k] )
           diff_f_sum += abs(kfile['profile_ext']['ffprim_ext'][k])
       print ('diff_p =', diff_p/diff_p_sum)
       print ('diff_f =', diff_f/diff_f_sum)
       if diff_p/diff_p_sum < max_error and diff_f/diff_f_sum < max_error: iconv = 1
    return iconv

def fixbdry_kfile_init(shot, time, f_inefit):
    #-------------------------------------------------------------------
    # read input
    #

    print (f_inefit)

    inefit = Namelist(f_inefit,"r")
    rho_in = inefit["inefit"]["rho"]

    ip = inefit["inefit"]['ip'][0]
    b0 = inefit["inefit"]['b0'][0]
    r0  = inefit["inefit"]['r0'][0]

    npsi = len(rho_in) #129
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

    return 0

def scale_kfile(shot, time, Rs, Bs):

    kfile = Namelist("k%06d.%05d"%(shot,time))
    # shutil.copyfile("k%06d.%05d"%(shot,time), "k%06d.%05d_s"%(shot,time))

    kfile['in1']['plasma'] = array(kfile['in1']['plasma'])/(Rs*Bs)
    kfile['in1']['btor'  ] = array(kfile['in1']['btor'  ])/Bs
    kfile['in1']['rcentr'] = array(kfile['in1']['rcentr'])/Rs
    kfile['in1']['rzero' ] = array(kfile['in1']['rzero' ])/Rs
    kfile['in1']['xlim'  ] = array(kfile['in1']['xlim'  ])/Rs
    kfile['in1']['ylim'  ] = array(kfile['in1']['ylim'  ])/Rs
    kfile['in1']['rbdry' ] = array(kfile['in1']['rbdry' ])/Rs
    kfile['in1']['zbdry' ] = array(kfile['in1']['zbdry' ])/Rs
    kfile['in1']['prbdry'] = array(kfile['in1']['prbdry'])/Bs**2

    kfile['profile_ext']['pprime_ext'] = array(kfile['profile_ext']['pprime_ext'])/Bs*Rs**2
    kfile['profile_ext']['ffprim_ext'] = array(kfile['profile_ext']['ffprim_ext'])/Bs

    kfile.write("k%06d.%05d_s"%(shot,time))

########################################################################
#   UTILS

def inverse3x3(x, y):
    a = array(
      [ [x[0]**2, x[0], 1.0],
        [x[1]**2, x[1], 1.0],
        [x[2]**2, x[2], 1.0] ] )
    a_inv = linalg.inv(a)
    return dot(a_inv,y)

def cubic(p, x0):
    return p[0]*x0**2+p[1]*x0+p[2]

def extrapolation(rho, f, x0, y0):
    spl = zinterp(rho,f,s=0)

    if y0 == None:
        #p = inverse3x3(x0, [spl(x) for x in x0])
        tmp = spl(x0, der=0)
        p = inverse3x3(x0, tmp)
        for k in range(len(rho)):
            if rho[k] < x0[0]:
                f[k] = cubic(p, rho[k])
    else:
        #p = inverse3x3(x0, [y0, spl(x0[1]), spl(x0[2])])
        tmp = spl([x0[1], x0[2]], der=0)
        p = inverse3x3(x0, [y0, tmp[0], tmp[1]])
        for k in range(len(rho)):
            if rho[k] < x0[1]:
                f[k] = cubic(p, rho[k])

def extrapolation_sep(x0, y0):
    p = inverse3x3(x0, y0)
    return cubic(p, 1.0)

########################################################################
#  STATE IO

def io_input_from_instate(f_instate, f_inefit="inefit", mode='kinetic'):

    instate = Namelist(f_instate)['instate']

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
    ni = (zimp**2*(ne-nbfast-2.0*(nalpha+nhe))-zimp*(zeff*ne-nbfast-4.0*(nalpha+nhe)))/tmp;
    nz = (ne*(zeff-1.0)-2.0*(nalpha+nhe))/tmp;

    if mode == 'kinetic':
       print('mode = kinetic')
       pmhd = 1.602e3*(ne*te +(ni+nz)*ti)+2.0/3.0*1.0e6*(wbeam+walp)
    elif mode == 'pmhd':
       print('mode = pmhd')
       pmhd = instate['pmhd']
       if pmhd is None:
          pmhd = instate['p_eq'] 
    else:
       print('mode = total')
       #pmhd = instate['pmhd']
       pmhd = instate['p_eq']

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
    inefit.write(f_inefit)

def io_input_from_state(f_ps, f_inbc, f_inefit="inefit", mode="kinetic", ismooth=0, betan_target=-1.):
    # read inbc
    inbc = Namelist(f_inbc)["inbc"]
    if not inbc: inbc = Namelist(f_inbc)["instate"]

    # read plasma state
    r0 = inbc["r0"][0]
    b0 = abs(inbc["b0"][0])
    ip = inbc['ip'][0]*1.0e6

    ps = plasmastate('ips',1)
    ps.read(f_ps)

    nrho  = len(ps["rho"])
    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]*1.0e-19
#   ni    = ps["ni"][:]*1.0e-19
    ni    = sum(ps["ns"][1:,:], axis=0)*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    ne    = ps.cell2node(ne)
    ni    = ps.cell2node(ni)
    te    = ps.cell2node(te)
    ti    = ps.cell2node(ti)

    density_beam = ps.dump_profile(rho, "rho_nbi", "nbeami", k=0)*1.e-19
    wbeam = ps.dump_profile(rho, "rho_nbi", "eperp_beami", k=0) \
        + ps.dump_profile(rho, "rho_nbi", "epll_beami", k=0)
    wbeam = density_beam*wbeam

    density_alpha = ps.dump_profile(rho, "rho_fus", "nfusi", k=0)*1.e-19
    walpha = ps.dump_profile(rho, "rho_fus", "eperp_fusi", k=0) \
        + ps.dump_profile(rho, "rho_fus", "epll_fusi", k=0)
    walpha = density_alpha*walpha

    if mode == 'kinetic':
       # print('mode = kinetic')
       pmhd = 1.602e3*(ne*te+ni*ti)+2.0/3.0*1.602e3*(wbeam+walpha)
    else:
       # print('mode = total')
       pmhd = ps["P_eq"][:]

    jpar = ps.dump_j_parallel(rho, "rho_eq", "curt", r0, b0, tot=True)

    #---------------------

    vol  = ps["vol" ][:]
    a0  = ps["rMinor_mean"][-1]
    wtot = 1.5*1.602e3*(ne*te+ni*ti)+1.602e3*(wbeam+walpha)
    wtot = wtot*1.0e6

    vol  = ps["vol" ][:]
    mu0 = 4.0*pi*1.0e-7
    wsum = sum([(vol[i+1]-vol[i])*wtot[i] for i in range(len(wtot)-1)])
    betan = wsum/vol[-1]/1.5
    betan *= 2.0*mu0/b0**2
    betan /= fabs(ip/(a0*b0))
    betan *= 1.0e2
    print ('betan =',betan)

    if betan_target > 0 and betan > betan_target:
        print("betan scale", betan_target/betan)
        pmhd = pmhd*betan_target/betan

    #---------------------

    # rlim = ps["rlim"][:]
    # zlim = ps["zlim"][:]
    # nlim = len(rlim)

    rlim = inbc["rlim"][:]
    zlim = inbc["zlim"][:]
    nlim = len(rlim)
    print('nlim =', nlim)
    print('rlim =', rlim)

    # namelist for efit

    inefit = Namelist()
    inefit["inefit"]["ip"   ] = [inbc["ip"][0]*1.0e6]
    inefit["inefit"]["r0"   ] = inbc["r0"]
    inefit["inefit"]["b0"   ] = inbc["b0"]
    inefit["inefit"]["nrho" ] = [nrho]
    inefit["inefit"]["rho"  ] = rho
    inefit["inefit"]["press"] = pmhd
    # inefit["inefit"]["jpar" ] = abs(jpar)
    if jpar[0] < 0:
        inefit["inefit"]["jpar" ] = -jpar
    else:
        inefit["inefit"]["jpar" ] = jpar

    if nlim:
        inefit["inefit"]["nlim" ] = [nlim]
        inefit["inefit"]["rlim" ] = rlim
        inefit["inefit"]["zlim" ] = zlim
    else:
        inefit["inefit"]["nlim" ] = inbc["nlim"]
        inefit["inefit"]["rlim" ] = inbc["rlim"]
        inefit["inefit"]["zlim" ] = inbc["zlim"]
    inefit["inefit"]["nbdry"] = inbc["nbdry"]
    inefit["inefit"]["rbdry"] = inbc["rbdry"]
    inefit["inefit"]["zbdry"] = inbc["zbdry"]
    inefit.write(f_inefit)

def io_input_init(f_instate):
    instate = Namelist(f_instate)['instate']

    nrho  = instate['nrho'][0]
    rho = arange(nrho)/(nrho-1.0)

    jpar  = zeros(nrho)
    pmhd  = zeros(nrho)

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

def io_input_init_keyargs(**keyargs):

    nrho  = keyargs['nrho']
    rho = arange(nrho)/(nrho-1.0)

    inefit = Namelist()
    inefit["inefit"]["ip"   ] = [keyargs["ip"]*1.0e6]
    inefit["inefit"]["r0"   ] = [keyargs["r0"]]
    inefit["inefit"]["b0"   ] = [keyargs["b0"]]
    inefit["inefit"]["nrho" ] = [nrho]
    inefit["inefit"]["rho"  ] = rho
    inefit["inefit"]["press"] = zeros(nrho)
    inefit["inefit"]["jpar" ] = zeros(nrho)
    inefit["inefit"]["nlim" ] = [keyargs["nlim" ]]
    inefit["inefit"]["rlim" ] = [keyargs["rlim" ]]
    inefit["inefit"]["zlim" ] = [keyargs["zlim" ]]
    inefit["inefit"]["nbdry"] = [keyargs["nbdry"]]
    inefit["inefit"]["rbdry"] = [keyargs["rbdry"]]
    inefit["inefit"]["zbdry"] = [keyargs["zbdry"]]
    inefit.write("inefit")