"""
model equilibrium, profile adjust
"""

from component import Component

from numpy import *
from Namelist import Namelist
from fastranutil import namelist_default
from zmodelprof import profile_hat, profile_spline, profile_pedestal
from zinterp import zinterp
import formula as fml 
from instate_model import expand_profile

def update_state(kiter, f_instate, nmax_iter=100, const=None):
    #------------------------------------------------------------------
    # entry
    print ("\n"+72*"=")
    print ('= model driver: adjust profiles ')

    iconv = 0

    #------------------------------------------------------------------
    # read plasma state
    instate = Namelist(f_instate)
    for key in instate["instate"].keys():
        if key.upper() not in [ 'TOKAMAK_ID', 'PRESSURE_MODEL', 'CURRENT_MODEL' ]:
            instate["instate"][key] = array(instate["instate"][key])
        else:
            print (key)

    nrho = instate["instate"]["nrho"][0]

    r0 = instate["instate"]["r0"][0]
    b0 = abs(instate["instate"]["b0"][0])
    ip = instate["instate"]['ip'][0]

    vol = instate["inmetric"]["vol"]
    area = instate["inmetric"]["area"]
    ipol = instate["inmetric"]["ipol"]
#   a0 = instate["inmetric"]["aminor"][0]
    a0 = instate["inmetric"]["rminor"][-1]
    rhob = instate["inmetric"]["rhob" ][0]
    volp = array(instate["inmetric"]["volp"])
    g22 = array(instate["inmetric"]["g22"])

    rho = instate["instate"]["rho"]
    ne = instate["instate"]["ne"]
    ni = instate["instate"]["ni"]
    nz = instate["instate"]["nz"]
    te = instate["instate"]["te"]
    ti = instate["instate"]["ti"]

    density_beam = array(instate["instate"]["density_beam"])
    wbeam = array(instate["instate"]["wbeam"])

    j_bs = instate["instate"]["j_bs"]

    qmhd = instate["inmetric"]["qmhd"]
    pmhd = instate["inmetric"]["pmhd"]

    n_ion = instate["instate"]["n_ion"][0]
    z_ion = instate["instate"]["z_ion"]
    a_ion = instate["instate"]["a_ion"]
    f_ion = instate["instate"]["f_ion"]
    n_imp = instate["instate"]["n_imp"][0]
    z_imp = instate["instate"]["z_imp"]
    a_imp = instate["instate"]["a_imp"]
    f_imp = instate["instate"]["f_imp"]

    jpar_axis = instate["instate"]["jpar_axis"][0]
    jpar_sep = instate["instate"]["jpar_sep"][0]
    jpar_alpha = instate["instate"]["jpar_alpha"][0]
    jpar_beta = instate["instate"]["jpar_beta"][0]
    zeff = instate["instate"]["zeff"][0]
    nbeam_axis = instate["instate"]["nbeam_axis"][0]
    nbeam_sep = instate["instate"]["nbeam_sep"][0]
    nbeam_alpha = instate["instate"]["nbeam_alpha"][0]
    nbeam_beta = instate["instate"]["nbeam_beta"][0]
    tbeami = instate["instate"]["tbeami"][0]
    ptot_alpha = instate["instate"]["ptot_alpha"][0]
    ptot_beta = instate["instate"]["ptot_beta"][0]
    xmid = instate["instate"]["xmid"][0]
    xwid = instate["instate"]["xwid"][0]
    te_xmid = instate["instate"]["te_xmid"][0]
    te_xwid = instate["instate"]["te_xwid"][0]
    te_axis = instate["instate"]["te_axis"][0]
    te_ped = instate["instate"]["te_ped"][0]
    te_sep = instate["instate"]["te_sep"][0]
    te_alpha = instate["instate"]["te_alpha"][0]
    te_beta = instate["instate"]["te_beta"][0]
    ti_xmid = instate["instate"]["ti_xmid"][0]
    ti_xwid = instate["instate"]["ti_xwid"][0]
    ti_axis = instate["instate"]["ti_axis"][0]
    ti_ped = instate["instate"]["ti_ped"][0]
    ti_sep = instate["instate"]["ti_sep"][0]
    ti_alpha = instate["instate"]["ti_alpha"][0]
    ti_beta = instate["instate"]["ti_beta"][0]
    tbeami = instate["instate"]["tbeami"][0]

    ne_fit = instate["instate"]["ne_fit"][0]
    te_fit = instate["instate"]["te_fit"][0]
    ti_fit = instate["instate"]["ti_fit"][0]

    print ("SMOOTH :", ne_fit, te_fit, ti_fit)

    iterate_beam = namelist_default(instate, "instate", "iterate_beam", [1])[0]

    betan_th_target = instate["instate"]["betan_th"][0]
    betan_beam_target = instate["instate"]["betan_beam"][0]
    betan_target = instate["instate"]["betan"][0]
    pressure_model =  instate["instate"]["pressure_model"][0].strip().lower()
    print ("pressure_model = ", pressure_model)
    current_model = instate["instate"]["current_model"][0].strip().lower()
    print ("current_model = ", current_model)

    cboot = namelist_default(instate, "instate", "cboot", [1.0])[0]

    rho_jbdry = namelist_default(instate, "instate", "rho_jbdry", [0.85])[0]

    rho_jhat = namelist_default(instate, "instate", "rho_jhat", [0.5])[0]
    wid_jhat = namelist_default(instate, "instate", "wid_jhat", [0.2])[0]

    q0 = namelist_default(instate, "instate", "q0", [1.2])[0]
    j_alpha = namelist_default(instate, "instate", "j_alpha", [1.5])[0]

    rho_jpeak = namelist_default(instate, "instate", "rho_jpeak", [0])[0]
    jaxis = namelist_default(instate, "instate", "jaxis", [0])[0]
    jpeak = namelist_default(instate, "instate", "jpeak", [0])[0]
    jaxis1 = namelist_default(instate, "instate", "jaxis_prime", [0])[0]
    jbdry1 = namelist_default(instate, "instate", "jbdry_prime", [0])[0]

    if current_model == "para":
        print ('q0 = ', q0)
        print ('j_alpha = ', j_alpha)
        print ('rho_jbdry = ', rho_jbdry)

    if current_model == "hat":
        print ('rho_jhat = ', rho_jhat)
        print ('wid_jhat = ', wid_jhat)
        print ('rho_jbdry = ', rho_jbdry)

    #rout = instate["afile"]["rout"][0]
    bout = b0 #r0*b0/rout

    #------------------------------------------------------------------
    # calculate betan
    wth  = 1.5*1.602e3*(ne*te+(ni+nz)*ti)
    betan_th = fml.betan(w=wth, vol=vol, ip=ip, b0=bout, a0=a0)
    betan_beam = fml.betan(w=wbeam*1.0e6, vol=vol, ip=ip, b0=bout, a0=a0)

    print ('betan_th   =', betan_th)
    print ('betan_beam =', betan_beam)

    #betan_th_target =  betan_th_target - betan_beam    #--------------

    #------------------------------------------------------------------
    # scale : beam pressure
    pbeam = 2.0/3.0*wbeam
    if betan_beam_target > 0.0:
        density_beam = (nbeam_axis-nbeam_sep)*(1.0-rho**nbeam_alpha)**nbeam_beta+nbeam_sep
        wbeam = 3./2.*1.602e3*density_beam*tbeami
        betan_beam = fml.betan(w=wbeam, vol=vol, ip=ip, b0=bout, a0=a0)

        scale_beam_density = betan_beam_target/betan_beam
        density_beam *=scale_beam_density
        wbeam *= scale_beam_density
        pbeam = 2.0/3.0*wbeam

        print ('scaled to betan_beam_target = ', betan_beam_target)
        print ('density_beam(0) = ', density_beam[0])

    #------------------------------------------------------------------
    # charge balance
    a = sum(f_ion*z_ion)
    b = sum(f_imp*z_imp)
    c = sum(f_ion*z_ion)
    d = sum(f_imp*z_imp*z_imp)

    zne_adj = ne
    zzne_adj = ne*zeff

    zne_adj = zne_adj - 1.0*density_beam
    zzne_adj = zzne_adj - 1.0**2*density_beam

    nion = (zne_adj*d-zzne_adj*b)/(a*d-b*c)
    nimp = (zzne_adj*a-zne_adj*c)/(a*d-b*c)

    density_ion = array([f_ion[k]*nion for k in range(n_ion)])
    density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

    density_th = array([sum(tmp) for tmp in density_ion.transpose()])
    density_th+= array([sum(tmp) for tmp in density_imp.transpose()])

    ni = array([sum(tmp) for tmp in density_ion.transpose()])
    nz = array([sum(tmp) for tmp in density_imp.transpose()])

    #------------------------------------------------------------------
    # scale : thermal pressure
    if betan_th_target > 0.0:
        tscale_min = 0.0
        tscale_max = 10.0

        for iter_scale in range(nmax_iter):
            tscale = 0.5*(tscale_min+tscale_max)

            te = profile_pedestal(nrho, te_xmid, te_xwid, te_ped, te_sep, tscale*te_axis, te_alpha, te_beta, ifit=te_fit)(rho)
            ti = profile_pedestal(nrho, ti_xmid, ti_xwid, ti_ped, ti_sep, tscale*ti_axis, ti_alpha, ti_beta, ifit=ti_fit)(rho)

            pth  = 1.602e3*(ne*te+(ni+nz)*ti)
            # betan_th = cal_betan(1.5*pth,vol,ip=ip,b0=bout,a0=a0)
            betan_th = fml.betan(w=1.5*pth,vol=vol,ip=ip,b0=bout,a0=a0)

            if abs(betan_th-betan_th_target) < 0.001: break

            if betan_th > betan_th_target:
               tscale_max = tscale
            else:
               tscale_min = tscale

            print ("scale to match betan_th : %6.3e, %6.3e, %6.3e %d"%(tscale,betan_th,betan_th_target,nmax_iter))
    else:
        pth  = 1.602e3*(ne*te+(ni+nz)*ti)
        tscale = 1.0

    pth_ped = pth[where( rho>=xmid-xwid/2)[0][0]]

    print ('Temperature scale: %5.3f %5.3f %5.3f %5.3f'%(te[0],ti[0],betan_th_target/betan_th,tscale))

    #------------------------------------------------------------------
    # scale : total pressure
    if pressure_model == "kinetic":
        p_axis_min = 0.0
        p_axis_max = 1000.0e3

        pmhd = zeros(nrho)

        for iter_scale in range(nmax_iter):
            p_axis = 0.5*(p_axis_min+p_axis_max)

            xped = xmid-xwid/2
            for k,xval in enumerate(rho):
                if 1.0-xval/xped > 0.0: pmhd[k] = p_axis*(1.0-(xval/xped)**ptot_alpha)**ptot_beta + pth_ped
                else: pmhd [k] = pth[k]

            # betan = cal_betan(1.5*pmhd,vol,ip=ip,b0=bout,a0=a0)
            betan = fml.betan(w=1.5*pmhd,vol=vol,ip=ip,b0=bout,a0=a0)

            if abs(betan-betan_target) < 0.01: break

            if betan > betan_target:
               p_axis_max = p_axis
            else:
               p_axis_min = p_axis

            print ("scale to match betan : %6.3e,%6.3e,%6.3f"%(p_axis*1.0e-3,betan,betan_target))

        print ('p_axis = ',p_axis)
        print ('betan = ',betan)

    elif pressure_model == "total":
        ptot_axis    = instate["instate"]["ptot_axis" ][0]
        ptot_ped     = instate["instate"]["ptot_ped"  ][0]
        ptot_sep     = instate["instate"]["ptot_sep"  ][0]
        ptot_alpha   = instate["instate"]["ptot_alpha"][0]
        ptot_beta    = instate["instate"]["ptot_beta" ][0]

        p_axis_min = 0.0
        p_axis_max = 100.0*ptot_axis

        pmhd = zeros(nrho)

        for iter_scale in range(nmax_iter):
            p_axis = 0.5*(p_axis_min+p_axis_max)

            pmhd = profile_pedestal(nrho, xmid, xwid, ptot_ped, ptot_sep, p_axis, ptot_alpha, ptot_beta)(rho)

            betan = fml.betan(w=1.5*pmhd, vol=vol, ip=ip, b0=bout, a0=a0)

            if abs(betan-betan_target) < 0.01: break

            if betan > betan_target:
               p_axis_max = p_axis
            else:
               p_axis_min = p_axis

            print ("scale to match betan : %6.3e,%6.3e"%(p_axis*1.0e-3,betan))

        print ('p_axis = ',p_axis)
        print ('betan = ',betan)

    elif pressure_model == "full":
        pmhd = pth + pbeam

    else:
        raise Exception ("check pressumre_model")

    #------------------------------------------------------------------
    # scale : current

    ibdry = where(rho>=rho_jbdry)[0][0]
    jbdry = j_bs[ibdry]

    def total_current(rho_jhat):
        jpeak_min = 0.1
        jpeak_max = 10.0

        j_hat = profile_hat(nrho,rho_jhat,wid_jhat)
        j_core = j_hat(rho)-j_hat[ibdry]

        j_tot = zeros(nrho)
        curt = zeros(nrho)

        for iter_scale in range(nmax_iter):

            jpeak = 0.5*(jpeak_min+jpeak_max)

            for i in range(nrho):
                if rho[i]<rho_jbdry: j_tot[i] = jpeak*j_core[i]+jbdry
                else: j_tot[i] = j_bs[i]

            for i in range(nrho-1):
                dV = vol[i+1]-vol[i]
                jparm = 0.5*(j_tot[i]+j_tot[i+1])
                ipolm = 0.5*(ipol[i]+ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

            if abs(curt[-1]-ip) < 1.0e-5: break

            if curt[-1] > ip:
               jpeak_max = jpeak
            else:
               jpeak_min = jpeak

        rhob = instate["inmetric"]["rhob"][0]
        g22  = array(instate["inmetric"]["g22"])
        volp = array(instate["inmetric"]["volp"])

        drho = rhob/(nrho-1.)
        jpar0 = 0.5*(j_tot[0]+j_tot[1])
        volp0 = 0.5*(volp[0]+volp[1])
        ipol0 = 0.5*(ipol[0]+ipol[1])
        dpsidrho = 0.2/g22[1]*(jpar0*volp0*drho/ipol0**2)
        q0_cal = 2.0*pi*b0*rho[1]*rhob/dpsidrho

        return q0_cal, j_tot

    if current_model == "hat":
        if q0>0.0:
            if kiter > 0: dq0 = 0.5*(q0 - instate["afile"]["qaxis"])
            else: dq0 = 0.0

            rho_jhat_min = namelist_default(instate, "instate", "rho_jhat_min", [0.1])[0]
            rho_jhat_max = namelist_default(instate, "instate", "rho_jhat_max", [0.6])[0]

            for iter_scale in range(nmax_iter):
                rho_jhat = 0.5*(rho_jhat_min+rho_jhat_max)
                q0_cal, j_tot = total_current(rho_jhat)
                if abs(q0_cal-q0) < 0.001: break
                print ("%6.3f %6.3f %6.3f %6.3f"%(q0_cal, q0, rho_jhat_min, rho_jhat_max))

               #if q0_cal > q0+dq0:
                if q0_cal > q0:
                    rho_jhat_max = rho_jhat
                else:
                    rho_jhat_min = rho_jhat

            print ("%6.3f %6.3f %6.3f %6.3f: dq0=%6.3f"%(q0_cal, q0, rho_jhat_min, rho_jhat_max,dq0))

        else:
            q0_cal, j_tot = total_current(rho_jhat)

        const["rho_jhat"] = rho_jhat

    elif current_model == "broad":
        jpeak_min = 0.1
        jpeak_max = 10.0

        j_tot = zeros(nrho)
        curt = zeros(nrho)

        for iter_scale in range(nmax_iter):

            x = 0.5*(jpeak_min+jpeak_max)

            x0 = [0.0, rho_jpeak, rho_jbdry]
            y0 = [jaxis/jpeak*x, x, jbdry]
            jaxis1 = x*(1. - jaxis/jpeak)/rho_jpeak
            yp = [jaxis1, 0.0, jbdry1]
            spl = profile_spline(x=x0, y=y0, yp=yp)

            j_tot = zeros(nrho)
            for i in range(nrho):
                if rho[i]<rho_jbdry: j_tot[i] = spl(rho[i])
                else: j_tot[i] = j_bs[i]

            curt = zeros(nrho)
            for i in range(nrho-1):
                dV = vol[i+1]-vol[i]
                jparm = 0.5*(j_tot[i]+j_tot[i+1])
                ipolm = 0.5*(ipol[i]+ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

            if abs(curt[-1]-ip) < 0.001: break

            if curt[-1] > ip:
               jpeak_max = x
            else:
               jpeak_min = x

        print('jpeak =', jpeak_min, jpeak_max)
        print('ip_cal = ', curt[-1], ip)

    elif current_model == "para":
        beta_min = 0.1
        beta_max = 5.0

        j_tot = zeros(nrho)
        curt = zeros(nrho)

        drho = rhob/(nrho-1.)
        volp0 = 0.5*(volp[0]+volp[1])
        ipol0 = 0.5*(ipol[0]+ipol[1])
        jpar0 = 2.0*pi*b0*rho[1]*rhob/q0
        jpar0 /= 0.2/g22[1]
        jpar0 /= volp0*drho/ipol0**2

        for iter_scale in range(100):

            beta = 0.5*(beta_min+beta_max)
            j_core = jpar0*(1.0-(rho/rho_jbdry)**j_alpha)**beta + jbdry

            for i in range(nrho):
                if rho[i]<rho_jbdry: j_tot[i] = j_core[i]
                else: j_tot[i] = j_bs[i]

            curt = zeros(nrho)
            for i in range(nrho-1):
                dV = vol[i+1]-vol[i]
                jparm = 0.5*(j_tot[i]+j_tot[i+1])
                ipolm = 0.5*(ipol[i]+ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

            if abs(curt[-1]-ip) < 0.001: break

            if curt[-1] > ip:
               beta_min = beta
            else:
               beta_max = beta

            print ('j0, beta, ip_cal, ip_target =',jpar0, beta, curt[-1], ip)

    else:
        raise Exception ("check current_model")

    instate["instate"]["pmhd"] = pmhd
    instate["instate"]["j_tot"] = j_tot
    instate["instate"]["te"] = te
    instate["instate"]["ti"] = ti
    if betan_beam_target > 0.0:
        instate["instate"]["density_beam"] =  density_beam
        instate["instate"]["wbeam"] = 1.0e-6*wbeam
    instate["instate"]["ni"] = ni
    instate["instate"]["nz"] = nz

    instate["instate"]["te_axis"]=[te[0]]
    instate["instate"]["ti_axis"]=[ti[0]]
    instate["instate"]["nbeam_axis"]=[density_beam[0]]

    instate.write(f_instate)

    iconv = 0
    return iconv, const

def constraint_pedestal_width(f_instate):
    #-- from efit

    instate = Namelist(f_instate)
    for key in instate["instate"].keys():
        if key.upper() not in [ 'TOKAMAK_ID', 'PRESSURE_MODEL', 'CURRENT_MODEL' ]:
            instate["instate"][key] = array(instate["instate"][key])

    pmhd_in = instate["instate"]["pmhd"]
    rho_in = instate["instate"]["rho"]

    pmhd = instate["inmetric"]["pmhd"]
    rho = instate["inmetric"]["rho"]
    psi = array(instate["inmetric"]["psi"])
    psi = psi/psi[-1]
    bp2 = instate["inmetric"]["bp2"]

    xmid = instate["instate"]["xmid"][0]
    xwid = instate["instate"]["xwid"][0]

    rho_ped = 1.0-xwid

    psi_ped = zinterp(rho,psi)(rho_ped)
    p_ped = zinterp(rho,pmhd)(rho_ped)

    p_ped_in = zinterp(rho_in,pmhd_in)(rho_ped)

    #-- pedestal betap
    mu0 = 4.*pi*1.e-7
    betap_ped = 2.0*mu0*p_ped/bp2[-1]
    xwid = 0.076*betap_ped**0.5
    xmid = 1.0-0.5*xwid - 0.01 #---------
    print ('adjust pedestal width:', xwid)

    instate["instate"]["xmid"] = [xmid]
    instate["instate"]["xwid"] = [xwid]
    instate["instate"]["ne_xmid"] = [xmid]
    instate["instate"]["ne_xwid"] = [xwid]
    instate["instate"]["te_xmid"] = [xmid]
    instate["instate"]["te_xwid"] = [xwid]
    instate["instate"]["ti_xmid"] = [xmid]
    instate["instate"]["ti_xwid"] = [xwid]

    expand_profile(instate)

    instate.write(f_instate)
