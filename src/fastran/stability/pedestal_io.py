"""
 -----------------------------------------------------------------------
 pedetal io
 -----------------------------------------------------------------------
"""

import numpy as np
from Namelist import Namelist
from fastran.util.modelprofile import profile_pedestal
from fastran.util.zinterp import zinterp
from fastran.util.loglinear import loglinear
from fastran.util.formula import get_ni
from fastran.util import formula
from fastran.util.fastranutil import namelist_default
from fastran.plasmastate.plasmastate import plasmastate


def write_input(fn_instate, fn_plasma_state, ps_backend):
    instate = Namelist(fn_instate)

    bt = namelist_default(instate, "instate", "b0", [-1.])[0]
    ip = namelist_default(instate, "instate", "ip", [-1.])[0]
    model_shape = namelist_default(instate, "instate", "model_shape", [0])[0]
    if model_shape in [0]:
        rb = instate["instate"]["rbdry"] 
        zb = instate["instate"]["zbdry"] 
        r = 0.5*( np.max(rb) + np.min(rb) )
        z = 0.5*( np.max(zb) + np.min(zb) )
        a = 0.5*( np.max(rb) - np.min(rb) )
        kappa = 0.5*( ( np.max(zb) - np.min(zb) )/ a )
        delta_u = ( r - rb[np.argmax(zb)] )/a
        delta_l = ( z - rb[np.argmin(zb)] )/a
        delta = 0.5*( delta_u + delta_l )
    elif model_shape in [1, 2]:
        r = namelist_default(instate, "instate", "r0", [-1.])[0] 
        a = namelist_default(instate, "instate", "a0", [-1.])[0] 
        kappa = namelist_default(instate, "instate", "kappa", [-1.])[0] 
        delta = namelist_default(instate, "instate", "delta", [-1.])[0]
    else:
        raise Exception("model_shape = {}, not in [0, 1, 2]".format(model_shape))
    neped = instate["instate"]["ne_ped"][0]
    nesep = instate["instate"]["ne_sep"][0]

    xwid = instate["instate"]["xwid"][0]
#   xmid = instate["instate"]["xmid"][0]
    xmid = 1. - 0.5*xwid
    rho_ped = xmid - 0.5*xwid

    #betan = instate["instate"]["betan"][0]
    if ps_backend == "INSTATE":
        rho = np.array(instate["instate"]["rho"])
        zeff = np.array(instate["instate"]["zeff"])
        zeff_ped = zinterp(rho, zeff)(rho_ped)

        ne = np.array(instate["instate"]["ne"])
        ni = np.array(instate["instate"]["ni"])
        nz = np.array(instate["instate"]["nz"])
        te = np.array(instate["instate"]["te"])
        ti = np.array(instate["instate"]["ti"])
        wbeam = np.array(instate["instate"]["wbeam"])
        vol = np.array(instate["inmetric"]["vol"])

        wth  = 1.5*1.602e3*(ne*te+(ni+nz)*ti)
        betan_th = formula.betan(w=wth, vol=vol, ip=ip, b0=bt, a0=a)
        betan_beam = formula.betan(w=wbeam*1.0e6, vol=vol, ip=ip, b0=bt, a0=a)
        betan = betan_th + betan_beam

    elif ps_backend == "PS":
        ps = plasmastate('ips', 1)
        ps.read(fn_plasma_state)
    
        rho = ps["rho"][:]
        zeff = ps.cell2node(ps["Zeff"][:])
        zeff_ped = zinterp(rho, zeff)(rho_ped)

        ne = ps.cell2node(ps["ns"][0, :]*1.0e-19)
        ni = ps.cell2node(np.sum(ps["ns"][1:, :], axis=0)*1.0e-19)
        te = ps.cell2node(ps["Ts"][0, :])
        ti = ps.cell2node(ps["Ti"][:])

        density_beam = ps.dump_profile(rho, "rho_nbi", "nbeami", k=0)*1.e-19
        wbeam = ps.dump_profile(rho, "rho_nbi", "eperp_beami", k=0) \
            + ps.dump_profile(rho, "rho_nbi", "epll_beami", k=0)
        wbeam = density_beam*wbeam
    
        density_alpha = ps.dump_profile(rho, "rho_fus", "nfusi", k=0)*1.e-19
        walpha = ps.dump_profile(rho, "rho_fus", "eperp_fusi", k=0) \
            + ps.dump_profile(rho, "rho_fus", "epll_fusi", k=0)
        walpha = density_alpha*walpha
    
        vol = zinterp(ps["rho_eq"], ps["vol"])(rho)
        wtot = 1.5*1.602e3*(ne*te + ni*ti) + 1.602e3*(wbeam + walpha)
        betan = formula.betan(w=wtot, vol=vol, ip=ip, b0=bt, a0=a)
       
    input_eped = Namelist()
    input_eped["input_eped"]["r"] = [r]
    input_eped["input_eped"]["a"] = [a]
    input_eped["input_eped"]["kappa"] = [kappa]
    input_eped["input_eped"]["delta"] = [delta]
    input_eped["input_eped"]["bt"] = [bt]
    input_eped["input_eped"]["ip"] = [ip] 
    input_eped["input_eped"]["neped"] = [neped] 
    input_eped["input_eped"]["nesep"] = [nesep] 
    input_eped["input_eped"]["zeffped"] = [zeff_ped] 
    input_eped["input_eped"]["betan"] = [betan]
    input_eped.write("input_eped")

def update_state(fn_inped, fn_instate, fn_plasma_state, ps_backend="INSTATE"):
    wped_model = loglinear(fn_inped, 'wped')
    pped_model = loglinear(fn_inped, 'pped')

    input_eped = Namelist("input_eped")

    wped = wped_model.get_(input_eped["input_eped"])
    pped = pped_model.get_(input_eped["input_eped"])

    ne_ped = input_eped["input_eped"]["neped"][0]
    zeff_ped = input_eped["input_eped"]["zeffped"][0]

    ni_ped, nz_ped = get_ni(ne_ped, zeff=zeff_ped)
    te_ped = 0.5*pped/(1.602*ne_ped)
    ti_ped = 0.5*pped/(1.602*(ni_ped+nz_ped))

    rho_ped = 1. - wped
    rho_top = 1. - 1.5*wped

    print('wped = {}'.format(wped))
    print('pped = {}'.format(pped))
    print('ne_ped = {}'.format(ne_ped))
    print('ni_ped = {}'.format(ni_ped))
    print('nz_ped = {}'.format(nz_ped))
    print('te_ped = {}'.format(te_ped))
    print('ti_ped = {}'.format(ti_ped))
    print('rho_ped = {}'.format(rho_ped))
    print('rho_top = {}'.format(rho_top))

    #-- input profiles from the previous iteration/timestep
    instate = Namelist(fn_instate)

    wped_in = instate["instate"]["xwid"][0]
    rho_ped_in = 1. - wped_in
    rho_top_in = 1. - 1.5*wped_in 

    y_in = {}
    if ps_backend == "INSTATE":
        rho = instate["instate"]["rho"]
        nrho = instate["instate"]["nrho"][0]
        for key in ["te", "ti", "ne"]:
            y_in[key] = instate["instate"][key]
    elif ps_backend == "PS":
        ps = plasmastate('ips', 1)
        ps.read(fn_plasma_state)
        rho   = ps["rho"][:]
        nrho = len(rho)
        y_in["ne"] = ps.cell2node_bdry(ps["ns"][0,:])*1.0e-19
        y_in["te"] = ps.cell2node_bdry(ps["Ts"][0,:])
        y_in["ti"] = ps.cell2node_bdry(ps["Ti"][:])

    profile_input = {}
    for key in ["te", "ti", "ne"]:
        profile_input[key] = zinterp(rho, y_in[key]) 
    
    ysep_in, ytop_in, yaxis_in = {}, {}, {}
    for key in ["te", "ti", "ne"]:
        ysep_in[key] = instate["instate"][key+"_sep"][0]
        ytop_in[key] = profile_input[key](rho_top_in)
        yaxis_in[key] = y_in[key][0]

    te_pedestal = profile_pedestal(nrho, 1.0-0.5*wped, wped, te_ped, ysep_in["te"], yaxis_in["te"], alpha=1.5, beta=1.5)
    ti_pedestal = profile_pedestal(nrho, 1.0-0.5*wped, wped, ti_ped, ysep_in["ti"], yaxis_in["ti"], alpha=1.5, beta=1.5)
    ne_pedestal = profile_pedestal(nrho, 1.0-0.5*wped, wped, ne_ped, ysep_in["ne"], yaxis_in["ne"], alpha=1.5, beta=1.5)

    te_top = te_pedestal(rho_top)
    ti_top = ti_pedestal(rho_top)
    ne_top = ne_pedestal(rho_top)

    te = np.zeros(nrho)
    ti = np.zeros(nrho)
    ne = np.zeros(nrho)
    for k in range(nrho):
        if rho[k] < rho_top:
            te[k] = profile_input["te"](rho[k] * rho_top_in/rho_top) + te_top - ytop_in["te"]
            ti[k] = profile_input["ti"](rho[k] * rho_top_in/rho_top) + ti_top - ytop_in["ti"]
            ne[k] = profile_input["ne"](rho[k] * rho_top_in/rho_top) + ne_top - ytop_in["ne"]
        else:
            te[k] = te_pedestal(rho[k])
            ti[k] = ti_pedestal(rho[k])
            ne[k] = ne_pedestal(rho[k])

    instate["instate"]["xwid"] = [wped]
    instate["instate"]["xmid"] = [1. - 0.5*wped]

    instate["instate"]["te"] = te
    instate["instate"]["te_ped"] = [te_ped]
    instate["instate"]["te_xwid"] = [wped] 
    instate["instate"]["te_xmid"] = [1. - 0.5*wped] 

    instate["instate"]["ti"] = ti
    instate["instate"]["ti_ped"] = [ti_ped]
    instate["instate"]["ti_xwid"] = [wped] 
    instate["instate"]["ti_xmid"] = [1. - 0.5*wped] 

    instate["instate"]["ne"] = ne
    instate["instate"]["ne_ped"] = [ne_ped]
    instate["instate"]["ne_xwid"] = [wped] 
    instate["instate"]["ne_xmid"] = [1. - 0.5*wped] 

    if ps_backend == "INSTATE":
        instate["instate"]["ne"] = ne
        instate["instate"]["te"] = te
        instate["instate"]["ti"] = ti
    elif ps_backend == "PS":
        ps["ns"][0] = 1.0e19*ps.node2cell(ne)
        ps["Ts"][0,:] = ps.node2cell(te)
        nspec_th = len(ps["Ts"]) - 1
        for k in range(nspec_th):
            ps["Ts"][k+1,:] = ps.node2cell(ti)
        ps["Ti"][:] = ps.node2cell(ti)
        ps.store(fn_plasma_state)

    instate.write(fn_instate)
