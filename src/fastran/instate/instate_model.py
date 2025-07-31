"""
 -----------------------------------------------------------------------
 fastran instate file utility
 -----------------------------------------------------------------------
"""
from numpy import *
from Namelist import Namelist
from fastran.util.fastranutil import namelist_default
from fastran.util.zinterp import zinterp
from fastran.util.model_profile import profile_pedestal
from fastran.util.formula import get_ni, mu0
from fastran.util.shape_io import boundaryShape, set_shape

def expand_profile(instate):
    for key in instate["instate"].keys(): instate["instate"][key] = array(instate["instate"][key])
    nrho       = instate["instate"]["nrho"       ][0]
    n_ion      = instate["instate"]["n_ion"      ][0]
    z_ion      = instate["instate"]["z_ion"      ]
    a_ion      = instate["instate"]["a_ion"      ]
    f_ion      = instate["instate"]["f_ion"      ]
    n_imp      = instate["instate"]["n_imp"      ][0]
    z_imp      = instate["instate"]["z_imp"      ]
    a_imp      = instate["instate"]["a_imp"      ]
    f_imp      = instate["instate"]["f_imp"      ]
    ne_axis    = instate["instate"]["ne_axis"    ][0]
    ne_ped     = instate["instate"]["ne_ped"     ][0]
    ne_sep     = instate["instate"]["ne_sep"     ][0]
    ne_alpha   = instate["instate"]["ne_alpha"   ][0]
    ne_beta    = instate["instate"]["ne_beta"    ][0]
    ne_xmid    = instate["instate"]["ne_xmid"    ][0]
    ne_xwid    = instate["instate"]["ne_xwid"    ][0]
    te_axis    = instate["instate"]["te_axis"    ][0]
    te_ped     = instate["instate"]["te_ped"     ][0]
    te_sep     = instate["instate"]["te_sep"     ][0]
    te_alpha   = instate["instate"]["te_alpha"   ][0]
    te_beta    = instate["instate"]["te_beta"    ][0]
    te_xmid    = instate["instate"]["te_xmid"    ][0]
    te_xwid    = instate["instate"]["te_xwid"    ][0]
    ti_axis    = instate["instate"]["ti_axis"    ][0]
    ti_ped     = instate["instate"]["ti_ped"     ][0]
    ti_sep     = instate["instate"]["ti_sep"     ][0]
    ti_alpha   = instate["instate"]["ti_alpha"   ][0]
    ti_beta    = instate["instate"]["ti_beta"    ][0]
    ti_xmid    = instate["instate"]["ti_xmid"    ][0]
    ti_xwid    = instate["instate"]["ti_xwid"    ][0]
    omega_axis = instate["instate"]["omega_axis" ][0]
    omega_sep  = instate["instate"]["omega_sep"  ][0]
    omega_alpha= instate["instate"]["omega_alpha"][0]
    omega_beta = instate["instate"]["omega_beta" ][0]
    zeff_axis  = instate["instate"]["zeff_axis"  ][0]
    nbeam_axis = instate["instate"]["nbeam_axis" ][0]
    nbeam_sep  = instate["instate"]["nbeam_sep"  ][0]
    nbeam_alpha= instate["instate"]["nbeam_alpha"][0]
    nbeam_beta = instate["instate"]["nbeam_beta" ][0]
    tbeami     = instate["instate"]["tbeami"     ][0]
    xmid       = instate["instate"]["xmid"       ][0]
    xwid       = instate["instate"]["xwid"       ][0]
    ne_fit     = instate["instate"]["ne_fit"     ][0]
    te_fit     = instate["instate"]["te_fit"     ][0]
    ti_fit     = instate["instate"]["ti_fit"     ][0]

    rho = arange(nrho)/(nrho-1.0)

    ne = profile_pedestal(nrho, ne_xmid, ne_xwid, ne_ped, ne_sep, ne_axis, ne_alpha, ne_beta, ifit=ne_fit)(rho)
    te = profile_pedestal(nrho, te_xmid, te_xwid, te_ped, te_sep, te_axis, te_alpha, te_beta, ifit=te_fit)(rho)
    ti = profile_pedestal(nrho, ti_xmid, ti_xwid, ti_ped, ti_sep, ti_axis, ti_alpha, ti_beta, ifit=ti_fit)(rho)

    omega = (omega_axis-omega_sep)*(1.0-rho**omega_alpha)**omega_beta+omega_sep
    zeff = zeff_axis*ones(nrho)
    density_beam = (nbeam_axis-nbeam_sep)*(1.0-rho**nbeam_alpha)**nbeam_beta+nbeam_sep
    density_alpha = zeros(nrho)

    a = sum(f_ion*z_ion)
    b = sum(f_imp*z_imp)
    c = sum(f_ion*z_ion)
    d = sum(f_imp*z_imp*z_imp)

    zne_adj = ne
    zzne_adj = ne*zeff

    zne_adj = zne_adj - 1.0*density_beam
    zzne_adj = zzne_adj - 1.0**2*density_beam

    nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
    nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

    density_ion = array([f_ion[k]*nion for k in range(n_ion)])
    density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

    density_th = array([sum(tmp) for tmp in density_ion.transpose()])
    density_th+= array([sum(tmp) for tmp in density_imp.transpose()])

    ni = array([sum(tmp) for tmp in density_ion.transpose()])
    nz = array([sum(tmp) for tmp in density_imp.transpose()])

    instate["instate"]["rho"] = rho
    instate["instate"]["ne"] = ne
    instate["instate"]["ni"] = ni
    instate["instate"]["nz"] = nz
    instate["instate"]["te"] = te
    instate["instate"]["ti"] = ti
    instate["instate"]["zeff"] = zeff
    instate["instate"]["omega"] = omega
    instate["instate"]["density_beam"] = density_beam
    instate["instate"]["wbeam"] = 3./2.*1.602e3*density_beam*tbeami*1.0e-6
    instate["instate"]["density_alpha"] = zeros(nrho)
    instate["instate"]["walpha"] = zeros(nrho)

def instate_model(f_instate):
    #--- read instate
    print("instate_model started")
    instate = Namelist(f_instate)
    for key in list(instate["instate"].keys()):
        if key.upper() not in [ 'TOKAMAK_ID', 'PRESSURE_MODEL', 'CURRENT_MODEL' ]:
            instate["instate"][key] = array(instate["instate"][key])
    # for key in instate["instate"].keys():
    #     instate["instate"][key] = array(instate["instate"][key])

    nrho       = instate["instate"]["nrho"       ][0]
    n_ion      = instate["instate"]["n_ion"      ][0]
    z_ion      = instate["instate"]["z_ion"      ]
    a_ion      = instate["instate"]["a_ion"      ]
    f_ion      = instate["instate"]["f_ion"      ]
    n_imp      = instate["instate"]["n_imp"      ][0]
    z_imp      = instate["instate"]["z_imp"      ]
    a_imp      = instate["instate"]["a_imp"      ]
    f_imp      = instate["instate"]["f_imp"      ]
    ne_axis    = instate["instate"]["ne_axis"    ][0]
    ne_ped     = instate["instate"]["ne_ped"     ][0]
    ne_sep     = instate["instate"]["ne_sep"     ][0]
    ne_alpha   = instate["instate"]["ne_alpha"   ][0]
    ne_beta    = instate["instate"]["ne_beta"    ][0]
    ne_xmid    = instate["instate"]["ne_xmid"    ][0]
    ne_xwid    = instate["instate"]["ne_xwid"    ][0]
    ne_fit     = instate["instate"]["ne_fit"     ][0]
    te_axis    = instate["instate"]["te_axis"    ][0]
    te_ped     = instate["instate"]["te_ped"     ][0]
    te_sep     = instate["instate"]["te_sep"     ][0]
    te_alpha   = instate["instate"]["te_alpha"   ][0]
    te_beta    = instate["instate"]["te_beta"    ][0]
    te_xmid    = instate["instate"]["te_xmid"    ][0]
    te_xwid    = instate["instate"]["te_xwid"    ][0]
    te_fit     = instate["instate"]["te_fit"     ][0]
    ti_axis    = instate["instate"]["ti_axis"    ][0]
    ti_ped     = instate["instate"]["ti_ped"     ][0]
    ti_sep     = instate["instate"]["ti_sep"     ][0]
    ti_alpha   = instate["instate"]["ti_alpha"   ][0]
    ti_beta    = instate["instate"]["ti_beta"    ][0]
    ti_xmid    = instate["instate"]["ti_xmid"    ][0]
    ti_xwid    = instate["instate"]["ti_xwid"    ][0]
    ti_fit     = instate["instate"]["ti_fit"     ][0]
    omega_axis = instate["instate"]["omega_axis" ][0]
    omega_sep  = instate["instate"]["omega_sep"  ][0]
    omega_alpha= instate["instate"]["omega_alpha"][0]
    omega_beta = instate["instate"]["omega_beta" ][0]
    jpar_axis  = instate["instate"]["jpar_axis"  ][0]
    jpar_sep   = instate["instate"]["jpar_sep"   ][0]
    jpar_alpha = instate["instate"]["jpar_alpha" ][0]
    jpar_beta  = instate["instate"]["jpar_beta"  ][0]
    zeff_axis  = instate["instate"]["zeff_axis"  ][0]
    nbeam_axis = instate["instate"]["nbeam_axis" ][0]
    nbeam_sep  = instate["instate"]["nbeam_sep"  ][0]
    nbeam_alpha= instate["instate"]["nbeam_alpha"][0]
    nbeam_beta = instate["instate"]["nbeam_beta" ][0]
    tbeami     = instate["instate"]["tbeami"     ][0]
    r0         = instate["instate"]["r0"         ][0]
    ip         = instate["instate"]["ip"         ][0]
    b0         = instate["instate"]["b0"         ][0]
    xmid       = instate["instate"]["xmid"       ][0]
    xwid       = instate["instate"]["xwid"       ][0]

    b0 = abs(b0)

    #--- boundary shape
    model_shape = namelist_default(instate,"instate","model_shape",[0])[0]
    print ('model_shape:', model_shape)

    if model_shape == 1:
        rb, zb, rlim, zlim = \
            set_shape(
                R0 = instate["instate"]["r0"][0],
                a0 = instate["instate"]["a0"][0],
                kappa = instate["instate"]["kappa"][0],
                delta = instate["instate"]["delta"][0],
                nt = instate["instate"]["nbdry"][0])
    elif model_shape == 2:
        print ("Luce Shape")
        R0 = instate["instate"]["r0"][0]
        a0 = instate["instate"]["a0"][0]
        eps = a0/R0
        kapu = instate["instate"]["kappa"][0]
        kapl = instate["instate"]["kappa"][0]
        delu = instate["instate"]["delta"][0]
        dell = instate["instate"]["delta"][0]
        z0 = namelist_default(instate,"instate","z0",[0])[0]
        np = instate["instate"]["nbdry"][0]

        zetaou = 0.
        zetaiu = 0.
        zetail = 0.
        zetaol = 0.
        zoffset = 0.

        rb, zb, zref = boundaryShape(a0, eps, kapu, kapl, delu, dell, zetaou, zetaiu, zetail, zetaol, zoffset,
                          upnull=True, lonull=True, npts=np, doPlot=False)

        rb = append(rb, rb[0])
        zb = append(zb, zb[0])

        zb += z0

        dlim = 0.05
        rmax = max(rb) + dlim
        rmin = min(rb) - dlim
        zmax = max(zb) + dlim
        zmin = min(zb) - dlim
        rlim = [ rmax, rmin, rmin, rmax, rmax ]
        zlim = [ zmax, zmax, zmin, zmin, zmax ]
    else:
        rb, zb, rlim, zlim = \
            instate["instate"]["rbdry"], \
            instate["instate"]["zbdry"], \
            instate["instate"]["rlim"], \
            instate["instate"]["zlim"]

    instate["instate"]["nbdry"] = [len(rb)]
    instate["instate"]["rbdry"] = rb
    instate["instate"]["zbdry"] = zb
    instate["instate"]["nlim"] = [len(rlim)]
    instate["instate"]["rlim"] = rlim
    instate["instate"]["zlim"] = zlim

    #--- pedestal
    if xmid > 0:
        ne_xmid = te_xmid = ti_xmid =  xmid
    if xwid > 0:
        ne_xwid = te_xwid = ti_xwid =  xwid

    betan_ped = instate["instate"]["betan_ped"][0]
    print ('betan_ped = ',betan_ped)
    if betan_ped > 0.0:
        ne_xwid = te_xwid = ti_xwid = xwid
        ne_xmid = te_xmid = ti_xmid = xmid

        rb = array(instate["instate"]["rbdry"])
        zb = array(instate["instate"]["zbdry"])
        a0 = 0.5*( max(rb) - min(rb) )

        c_betan = 4.0*1.602e5*ne_ped*mu0/b0**2*(a0*b0/ip)
        te_ped = betan_ped/c_betan
        ti_ped = te_ped
        print ("PEDEDSTAL BETAN = ",betan_ped, te_ped)

        instate["instate"]["te_ped" ] = [te_ped ]
        instate["instate"]["te_xwid"] = [te_xwid]
        instate["instate"]["te_xmid"] = [te_xmid]

        instate["instate"]["ti_ped" ] = [ti_ped ]
        instate["instate"]["ti_xwid"] = [ti_xwid]
        instate["instate"]["ti_xmid"] = [ti_xmid]

    use_inped = instate["instate"]["use_inped"][0]
    if  use_inped == 1:
        print("use inped")
        nesep = ne_ped * instate["inped"]["nesep"][0]
        tesep_inped = instate["inped"]["tesep"][0]
        if tesep_inped > 0:
            tesep = tesep_inped
        else:
            tesep = abs(tesep_model)*teped
        tisep_inped = instate["inped"]["tisep"][0]
        if tisep_inped > 0:
            tisep = tisep_inped
        else:
            tisep = abs(tisep_inped)*ti_ped
        neped = ne_ped

        teped = instate["inped"]["teped_const"][0]
        teped*= ip**instate["inped"]["teped_ip"][0]
        teped*= b0**instate["inped"]["teped_bt"][0]
        teped*= neped**instate["inped"]["teped_neped"][0]
        teped*= nesep**instate["inped"]["teped_nesep"][0]

        wped = instate["inped"]["wped_const"][0]
        wped*= ip**instate["inped"]["wped_ip"][0]
        wped*= b0**instate["inped"]["wped_bt"][0]
        wped*= neped**instate["inped"]["wped_neped"][0]
        wped*= nesep**instate["inped"]["wped_nesep"][0]

        xmid = 1.0-0.5*wped
        xwid = wped

        niped, nzped = get_ni(neped, zeff=zeff_axis)
        tiped = teped * neped / (niped+nzped)
        ti0 = te_axis*tiped/teped

        print ('teped, tiped = ', teped, tiped)
        print ('niped, nzped = ', niped, nzped)
        print ('nesep = ', nesep)
        print ('tesep, tisep = ', tesep, tisep)

        ne_xmid = te_xmid = ti_xmid =  xmid
        ne_xwid = te_xwid = ti_xwid =  xwid

        ne_ped = neped
        ne_sep = nesep
        te_ped = teped
        ti_ped = teped

    print(">>>>> check te_ped = ", te_ped)
    print(">>>>> check te_xmid = ", te_xmid)
    print(">>>>> check te_xwid = ", te_xwid)

    #--- construnct profile
    rho = arange(nrho)/(nrho-1.0)
    ne = profile_pedestal(nrho, ne_xmid, ne_xwid, ne_ped, ne_sep, ne_axis, ne_alpha, ne_beta, ifit=ne_fit)(rho)
    te = profile_pedestal(nrho, te_xmid, te_xwid, te_ped, te_sep, te_axis, te_alpha, te_beta, ifit=te_fit)(rho)
    ti = profile_pedestal(nrho, ti_xmid, ti_xwid, ti_ped, ti_sep, ti_axis, ti_alpha, ti_beta, ifit=ti_fit)(rho)

    omega = (omega_axis-omega_sep)*(1.0-rho**omega_alpha)**omega_beta+omega_sep
    zeff = zeff_axis*ones(nrho)
    j_tot = (jpar_axis-jpar_sep)*(1.0-rho**jpar_alpha)**jpar_beta+jpar_sep
    pmhd = 1.0e3*(1.0-rho**1.5)**1.5
    density_beam = (nbeam_axis-nbeam_sep)*(1.0-rho**nbeam_alpha)**nbeam_beta+nbeam_sep
    density_alpha = zeros(nrho)
    tbeami = tbeami*ones(nrho)

    #--- density
    a = sum(f_ion*z_ion)
    b = sum(f_imp*z_imp)
    c = sum(f_ion*z_ion)
    d = sum(f_imp*z_imp*z_imp)

    zne_adj = ne
    zzne_adj = ne*zeff

    zne_adj = zne_adj - 1.0*density_beam
    zzne_adj = zzne_adj - 1.0**2*density_beam

    nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
    nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

    density_ion = array([f_ion[k]*nion for k in range(n_ion)])
    density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

    density_th = array([sum(tmp) for tmp in density_ion.transpose()])
    density_th+= array([sum(tmp) for tmp in density_imp.transpose()])

    ni = array([sum(tmp) for tmp in density_ion.transpose()])
    nz = array([sum(tmp) for tmp in density_imp.transpose()])

    #--- put to instate
    instate["instate"]["rho"] = rho
    instate["instate"]["ne"] = ne
    instate["instate"]["ni"] = ni
    instate["instate"]["nz"] = nz
    instate["instate"]["te"] = te
    instate["instate"]["ti"] = ti
    instate["instate"]["zeff"] = zeff
    instate["instate"]["omega"] = omega
    instate["instate"]["density_beam"] = density_beam
    instate["instate"]["wbeam"] = 3./2.*1.602e3*density_beam*tbeami*1.0e-6
    instate["instate"]["density_alpha"] = zeros(nrho)
    instate["instate"]["walpha"] = zeros(nrho)
    instate["instate"]["j_tot"] = j_tot
    instate["instate"]["pmhd"] = pmhd

    for k in range(n_ion):
        instate["check"]["ni_{}".format(k)] = density_ion[k]
    for k in range(n_imp):
        instate["check"]["nz_{}".format(k)] = density_imp[k]

    #--- zeros
    for key in [
        "j_oh", "j_bs", "j_nb", "j_ec", "j_ic", \
        "pe_nb", "pe_ec", "pe_ic", "pe_fus", "pe_ionization", "p_rad", \
        "pi_nb", "pi_ec", "pi_ic", "pi_fus", "pi_ionization", "pi_cx", "p_ohm", "p_ei", \
        "torque_nb", "torque_in", "se_nb", "se_ionization", "si_nb", "si_ionization", \
        "q", "psipol",  \
        "chie", "chii", "p_eq" ]:
        instate["instate"][key] = zeros(nrho)

    #--- write
    instate.write(f_instate)
