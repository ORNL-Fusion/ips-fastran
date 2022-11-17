"""
model equilibrium, profile adjust
"""

import numpy as np
from Namelist import Namelist
from fastran.util.modelprofile import profile_pedestal, profile_hat, profile_spline
from fastran.util.zinterp import zinterp
from fastran.instate.instate_model import expand_profile
from fastran.state.instate import Instate
import fastran.util.formula as fml
from ipsframework import Component


def update_state(kiter, f_instate, nmax_iter=100, const=None):
    instate = Instate(f_instate)

    r0 = instate["r0"][0]
    a0 = instate["a0"][0]
    b0 = abs(instate["b0"][0])
    ip = instate["ip"][0]
    bout = b0
    vol = instate.data["inmetric"]["vol"]
    rho = instate["rho"]
    nrho = instate["nrho"][0]

    # -------------------------------------------------------------------
    # calculate betan
    wth = 1.5*1.602e3*(instate["ne"]*instate["te"] + (instate["ni"] + instate["nz"])*instate["ti"])
    betan_th = fml.betan(w=wth, vol=vol, ip=ip, b0=bout, a0=a0)
    betan_beam = fml.betan(w=instate["wbeam"]*1.e6, vol=vol, ip=ip, b0=bout, a0=a0)

    print('betan_th   =', betan_th)
    print('betan_beam =', betan_beam)

    # -------------------------------------------------------------------
    # betan or h98 constraint
    x = instate.data["inmetric"]["rminor"]
    nebar = 0.
    for i in range(nrho - 1):
        nebar += 0.5*(instate["ne"][i + 1] + instate["ne"][i]) * (x[i+1] - x[i])
    nebar /= x[-1]
    print('nebar = ', nebar)

    kappa = instate["kappa"][0]
    m = np.sum( [ instate["a_ion"][k] * instate["f_ion"][k] for k in range(instate["n_ion"][0]) ] )
    print("m = ", m)

    vol_integrate = lambda vol, y: np.sum([(vol[i + 1] - vol[i]) * 0.5*(y[i + 1] + y[i]) for i in range(len(y) - 1)])

    pe_ec  = vol_integrate(vol, instate["pe_ec"]) #np.sum([(vol[i + 1] - vol[i])*instate["pe_ec"][i] for i in range(nrho - 1)])
    pe_ic  = vol_integrate(vol, instate["pe_ic"]) #np.sum([(vol[i + 1] - vol[i])*instate["pe_ic"][i] for i in range(nrho - 1)])
    pi_ic  = vol_integrate(vol, instate["pi_ic"]) #np.sum([(vol[i + 1] - vol[i])*instate["pi_ic"][i] for i in range(nrho - 1)])
    pe_nb  = vol_integrate(vol, instate["pe_nb"]) #np.sum([(vol[i + 1] - vol[i])*instate["pe_nb"][i] for i in range(nrho - 1)])
    pi_nb  = vol_integrate(vol, instate["pi_nb"]) #np.sum([(vol[i + 1] - vol[i])*instate["pi_nb"][i] for i in range(nrho - 1)])
    pe_fus = vol_integrate(vol, instate["pe_fus"]) #np.sum([(vol[i + 1] - vol[i])*instate["pe_fus"][i] for i in range(nrho - 1)])
    pi_fus = vol_integrate(vol, instate["pi_fus"]) #np.sum([(vol[i + 1] - vol[i])*instate["pi_fus"][i] for i in range(nrho - 1)])

    pinj = pe_ec + pe_ic + pi_ic + pe_nb + pi_nb + pe_fus + pi_fus

    print('pe_ec, pe_ic, pi_ic, pe_nb, pi_nb =', pe_ec, pe_ic, pi_ic, pe_nb, pi_nb)
    print('pe_fus, pi_fus =', pe_fus, pi_fus)
    print('pinj =', pinj)

    tau98 = 0.0562 \
          * ip**0.93 \
          * b0**0.15 \
          * pinj**-0.69 \
          * nebar**0.41 \
          * m**0.19 \
          * r0**1.97 \
          * (r0/a0)**-0.58 \
          * kappa**0.78

    print('tau98 =', tau98)

    h98_target = instate.default("h98_target", [-1.])[0]
    print('h98_target =', h98_target)
    if h98_target > 0:
        wth_target = tau98 * h98_target * pinj 
        mu0 = 4.*np.pi*1.e-7
        betan_th_target = 4./3.*mu0*a0/(ip*b0) * wth_target*1.0e6 / vol[-1] * 100.
        print('betan_th_target = ', betan_th_target) 
        betan_beam_target = instate.default("betan_beam_target", [0.])[0]
        betan_target = betan_th_target + betan_beam_target
        instate["betan_target"] = [betan_target]
    else:
        betan_target = instate["betan_target"][0]
        betan_beam_target = instate.default("betan_beam_target", [0.])[0]
        betan_th_target = betan_target - betan_beam_target

    # -------------------------------------------------------------------
    # scale : beam pressure

    if betan_beam_target > 0.:
        density_beam = (instate["nbeam_axis"][0] - instate["nbeam_sep"][0]) \
            * (1. - instate["rho"]**instate["nbeam_alpha"][0])**instate["nbeam_beta"][0] \
            + instate["nbeam_sep"][0]
        wbeam = 3./2.*1.602e3*density_beam*instate["tbeami"][0]
        betan_beam = fml.betan(w=wbeam, vol=vol, ip=ip, b0=bout, a0=a0)

        if betan_beam > 0:
            scale_beam_density = betan_beam_target/betan_beam
        else:
            raise Exception("*** beam ion density = 0")
        density_beam *= scale_beam_density
        wbeam *= scale_beam_density

        print('scaled to betan_beam_target = ', betan_beam_target)
        print('density_beam(0) = ', density_beam[0])

        instate["density_beam"] = density_beam
        instate["wbeam"] = wbeam
        instate["nbeam_axis"] = [density_beam[0]]

    # -------------------------------------------------------------------
    # particle and charge balance
    instate.particle_balance()

    # -------------------------------------------------------------------
    # scale : thermal pressure
    if betan_th_target > 0.:
        tscale_min = 0.
        tscale_max = 10.

        for iter_scale in range(nmax_iter):
            tscale = 0.5*(tscale_min + tscale_max)

            te = profile_pedestal(
                nrho,
                instate["te_xmid"][0], instate["te_xwid"][0], instate["te_ped"][0], instate["te_sep"][0],
                tscale*instate["te_axis"][0], instate["te_alpha"][0], instate["te_beta"][0],
                ifit=instate["te_fit"][0])(rho)
            ti = profile_pedestal(
                nrho,
                instate["ti_xmid"][0], instate["ti_xwid"][0], instate["ti_ped"][0], instate["ti_sep"][0],
                tscale*instate["ti_axis"][0], instate["ti_alpha"][0], instate["ti_beta"][0],
                ifit=instate["ti_fit"][0])(rho)

            pth = 1.602e3*(instate["ne"]*te + (instate["ni"] + instate["nz"])*ti)
            betan_th = fml.betan(w=1.5*pth, vol=vol, ip=ip, b0=bout, a0=a0)

            if abs(betan_th - betan_th_target) < 0.001:
                break

            if betan_th > betan_th_target:
                tscale_max = tscale
            else:
                tscale_min = tscale

            print("scale to match betan_th : %6.3e, %6.3e, %6.3e" % (tscale, betan_th, betan_th_target))

        instate["te"] = te
        instate["ti"] = ti

        instate["te_axis"] = [te[0]]
        instate["ti_axis"] = [ti[0]]

    else:
        pth = 1.602e3*(instate["ne"]*insate["te"] + (instate["ni"] + instate["nz"])*instate["ti"])
        tscale = 1.0

    print('Temperature scale: %5.3f %5.3f %5.3f %5.3f' % (te[0], ti[0], betan_th_target/betan_th, tscale))

    # ------------------------------------------------------------------
    # scale : total pressure
    pressure_model = instate.default("pressure_model", ["kinetic"])[0].lower()

    if pressure_model in ["kinetic", "full"]:
        pmhd = pth + 2./3.*instate["wbeam"]

    elif pressure_model in ["total"]:
        p_axis_min = 0.0
        p_axis_max = 100.0*instate["ptot_axis"][0]

        pmhd = np.zeros(nrho)

        for iter_scale in range(nmax_iter):
            p_axis = 0.5*(p_axis_min + p_axis_max)

            pmhd = profile_pedestal(
                nrho,
                instate["ptot_xmid"][0], instate["ptot_xwid"][0], instate["ptot_ped"][0], instate["ptot_sep"][0],
                p_axis, instate["ptot_alpha"][0], instate["ptot_beta"][0])(rho)

            betan = fml.betan(w=1.5*pmhd, vol=vol, ip=ip, b0=bout, a0=a0)

            if abs(betan - betan_target) < 0.01:
                break

            if betan > betan_target:
                p_axis_max = p_axis
            else:
                p_axis_min = p_axis

            print("scale to match betan : %6.3e, %6.3e" % (p_axis*1.0e-3, betan))

        print('p_axis = ', p_axis)
        print('betan = ', betan)

    else:
        raise Exception("check pressumre_model")

    instate["pmhd"] = pmhd

    # ------------------------------------------------------------------
    # current constraint
    current_model = instate.default("current_model", ["broad"])[0].lower()

    rho_jbdry = instate.default("rho_jbdry", [0.85])[0]
    ibdry = np.where(rho >= rho_jbdry)[0][0]
    jbdry = instate["j_bs"][ibdry]

    j_tot = np.zeros(nrho)
    curt = np.zeros(nrho)
    ipol = instate.data["inmetric"]["ipol"]
    volp = instate.data["inmetric"]["volp"]
    g22 = instate.data["inmetric"]["g22"]
    rhob = instate.data["inmetric"]["rhob"][0]

    if current_model == "broad":
        rho_jpeak = instate.default("rho_jpeak", [0.6])[0]
        jaxis = instate.default("jaxis", [0.7])[0]
        jaxis1 = instate.default("jaxis_prime", [1.])[0]
        jbdry1 = instate.default("jbdry_prime", [-1.])[0]

        jpeak_min, jpeak_max = 0.1, 10.0
        for iter_scale in range(nmax_iter):
            x = 0.5*(jpeak_min + jpeak_max)

            x0 = [0.0, rho_jpeak, rho_jbdry]
            y0 = [jaxis*x, x, jbdry]
            jaxis1 = x*(1. - jaxis)/rho_jpeak
            yp = [jaxis1, 0., jbdry1]
            spl = profile_spline(x=x0, y=y0, yp=yp)

            for i in range(nrho):
                if rho[i] < rho_jbdry:
                    j_tot[i] = spl(rho[i])
                else:
                    j_tot[i] = instate["j_bs"][i]

            for i in range(nrho - 1):
                dV = vol[i+1] - vol[i]
                jparm = 0.5*(j_tot[i] + j_tot[i+1])
                ipolm = 0.5*(ipol[i] + ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i] + jparm*dV/(2.0*np.pi*r0*ipolm**2))*ipol[i + 1]

            if abs(curt[-1] - ip) < 0.001:
                break

            if curt[-1] > ip:
                jpeak_max = x
            else:
                jpeak_min = x

        print('jpeak =', jpeak_min, jpeak_max)
        print('ip_cal = ', curt[-1], ip)

    elif current_model in ["parabolic", "para"]:
        q0 = instate.default("q0", [1.2])[0]
        j_alpha = instate.default("j_alpha", [1.5])[0]

        drho = rhob/(nrho - 1.)
        volp0 = 0.5*(volp[0] + volp[1])
        ipol0 = 0.5*(ipol[0] + ipol[1])
        jpar0 = 2.0*np.pi*b0*rho[1]*rhob/q0
        jpar0 /= 0.2/g22[1]
        jpar0 /= volp0*drho/ipol0**2

        j_beta_min, j_beta_max = 0.1, 5.0
        for iter_scale in range(nmax_iter):
            j_beta = 0.5*(j_beta_min + j_beta_max)
            j_core = jpar0*(1.0-(rho/rho_jbdry)**j_alpha)**j_beta + jbdry

            for i in range(nrho):
                if rho[i] < rho_jbdry:
                    j_tot[i] = j_core[i]
                else:
                    j_tot[i] = instate["j_bs"][i]

            for i in range(nrho - 1):
                dV = vol[i+1] - vol[i]
                jparm = 0.5*(j_tot[i] + j_tot[i+1])
                ipolm = 0.5*(ipol[i] + ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i] + jparm*dV/(2.0*np.pi*r0*ipolm**2))*ipol[i + 1]

            if abs(curt[-1] - ip) < 0.001:
                break

            if curt[-1] > ip:
                j_beta_min = j_beta
            else:
                j_beta_max = j_beta

            print('j0, j_beta, ip_cal, ip_target =', jpar0, j_beta, curt[-1], ip)

    elif current_model == "hat":
        rho_jhat = instate.default("rho_jhat", [0.5])[0]
        wid_jhat = instate.default("wid_jhat", [0.2])[0]
        rho_jhat_min = instate.default("rho_jhat_min", [0.1])[0]
        rho_jhat_max = instate.default("rho_jhat_max", [0.6])[0]
        q0 = instate.default("q0", [1.2])[0]

        for iter_scale_out in range(nmax_iter):
            rho_jhat = 0.5*(rho_jhat_min + rho_jhat_max)

            j_hat = profile_hat(nrho, rho_jhat, wid_jhat)
            j_core = j_hat(rho) - j_hat[ibdry]

            jpeak_min, jpeak_max = 0.1, 10.
            for iter_scale in range(nmax_iter):

                jpeak = 0.5*(jpeak_min + jpeak_max)

                for i in range(nrho):
                    if rho[i] < rho_jbdry:
                        j_tot[i] = jpeak*j_core[i]+jbdry
                    else:
                        j_tot[i] = instate["j_bs"][i]

                for i in range(nrho - 1):
                    dV = vol[i+1] - vol[i]
                    jparm = 0.5*(j_tot[i] + j_tot[i+1])
                    ipolm = 0.5*(ipol[i] + ipol[i+1])
                    curt[i+1] = (curt[i]/ipol[i] + jparm*dV/(2.0*np.pi*r0*ipolm**2))*ipol[i + 1]

                if abs(curt[-1] - ip) < 1.0e-5:
                    break

                if curt[-1] > ip:
                    jpeak_max = jpeak
                else:
                    jpeak_min = jpeak

            drho = rhob/(nrho - 1.)
            jpar0 = 0.5*(j_tot[0] + j_tot[1])
            volp0 = 0.5*(volp[0] + volp[1])
            ipol0 = 0.5*(ipol[0] + ipol[1])
            dpsidrho = 0.2/g22[1]*(jpar0*volp0*drho/ipol0**2)
            q0_cal = 2.0*np.pi*b0*rho[1]*rhob/dpsidrho

            if abs(q0_cal - q0) < 0.001:
                break
            print("%6.3f %6.3f %6.3f %6.3f" % (q0_cal, q0, rho_jhat_min, rho_jhat_max))

            if q0_cal > q0:
                rho_jhat_max = rho_jhat
            else:
                rho_jhat_min = rho_jhat

        print("%6.3f %6.3f %6.3f %6.3f" % (q0_cal, q0, rho_jhat_min, rho_jhat_max))

    else:
        raise Exception("check current_model")

    instate["j_tot"] = j_tot

    instate.write(f_instate)

    iconv = 0
    return iconv, const


def constraint_pedestal_width(f_instate):
    # -- from efit

    instate = Namelist(f_instate)
    for key in instate["instate"].keys():
        if key.upper() not in ['TOKAMAK_ID', 'PRESSURE_MODEL', 'CURRENT_MODEL']:
            instate["instate"][key] = np.array(instate["instate"][key])

    pmhd_in = instate["instate"]["pmhd"]
    rho_in = instate["instate"]["rho"]

    pmhd = instate["inmetric"]["pmhd"]
    rho = instate["inmetric"]["rho"]
    psi = np.array(instate["inmetric"]["psi"])
    psi = psi/psi[-1]
    bp2 = instate["inmetric"]["bp2"]

    xmid = instate["instate"]["xmid"][0]
    xwid = instate["instate"]["xwid"][0]

    rho_ped = 1.0-xwid

    psi_ped = zinterp(rho, psi)(rho_ped)
    p_ped = zinterp(rho, pmhd)(rho_ped)

    p_ped_in = zinterp(rho_in, pmhd_in)(rho_ped)

    # -- pedestal betap
    mu0 = 4.*np.pi*1.e-7
    betap_ped = 2.0*mu0*p_ped/bp2[-1]
    xwid = 0.076*betap_ped**0.5
    xmid = 1.0-0.5*xwid - 0.01  # ---------
    print('adjust pedestal width:', xwid)

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
