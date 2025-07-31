"""
 -----------------------------------------------------------------------
 pedetal io
 -----------------------------------------------------------------------
"""
import numpy as np
from scipy import optimize
from Namelist import Namelist
from fastran.util.model_profile import mtanh0, profile_pedestal
from fastran.util.zinterp import zinterp
from fastran.util.loglinear import loglinear
from fastran.util.formula import get_ni
from fastran.util import formula
from fastran.util.fastranutil import namelist_default
from fastran.plasmastate.plasmastate import plasmastate


class profile_core_pedestal:
    def __init__(self, x, y, wped, yped, yaxis, ysep, kind='mtanh'):
        self.x = x
        self.y = y
        self.spline = zinterp(x, y)

        if kind == 'mtanh':
            def f_yaxis(alpha):
                return mtanh0([1. - 0.5 * wped, wped, yped, ysep, alpha], 0.) - yaxis
            alpha0 = optimize.bisect(f_yaxis, 0., 1., xtol=1.0e-3)
            self.pedestal = zinterp(x, mtanh0([1. - 0.5 * wped, wped, yped, ysep, alpha0], x))
        elif kid == 'eped':
            self.pedestal = profile_pedestal(len(x), 1. - 0.5 * wped, wped, yped, ysep, yaxis, alpha=1.5, beta=1.5, ifit=0)

        self.wped = wped
        self.xtop = 1. - 1.5 * wped
        self.ytop = self.pedestal(self.xtop)

    def patch_eped(self, x):
        profile = np.zeros(len(x))
        for k in range(len(x)):
             if x[k] < self.xtop:
               # profile[k] = self.spline(x[k] * xtop/self.xtop) + self.ytop - self.spline(xtop)
                 profile[k] = self.spline(x[k]) + self.ytop - self.spline(self.xtop)
             else:
                 profile[k] = self.pedestal(x[k])
        return profile


def write_input(fn_instate, fn_plasma_state, ps_backend):
    instate = Namelist(fn_instate)

    bt = namelist_default(instate, 'instate', 'b0', [-1.])[0]
    ip = namelist_default(instate, 'instate', 'ip', [-1.])[0]
    model_shape = namelist_default(instate, 'instate', 'model_shape', [0])[0]
    if model_shape in [0]:
        rb = instate['instate']['rbdry']
        zb = instate['instate']['zbdry']
        r = 0.5*( np.max(rb) + np.min(rb) )
        z = 0.5*( np.max(zb) + np.min(zb) )
        a = 0.5*( np.max(rb) - np.min(rb) )
        kappa = 0.5*( ( np.max(zb) - np.min(zb) )/ a )
        delta_u = ( r - rb[np.argmax(zb)] )/a
        delta_l = ( z - rb[np.argmin(zb)] )/a
        delta = 0.5*( delta_u + delta_l )
    elif model_shape in [1, 2]:
        r = namelist_default(instate, 'instate', 'r0', [-1.])[0]
        a = namelist_default(instate, 'instate', 'a0', [-1.])[0]
        kappa = namelist_default(instate, 'instate', 'kappa', [-1.])[0]
        delta = namelist_default(instate, 'instate', 'delta', [-1.])[0]
    else:
        raise Exception('model_shape = {}, not in [0, 1, 2]'.format(model_shape))
    neped = instate['instate']['ne_ped'][0]
    nesep = instate['instate']['ne_sep'][0]
    zeff_ped = instate['instate']['zeff_ped'][0]

    # betan = instate['instate']['betan'][0]
    if ps_backend == 'INSTATE':
        rho = np.array(instate['instate']['rho'])
        ne = np.array(instate['instate']['ne'])
        ni = np.array(instate['instate']['ni'])
        nz = np.array(instate['instate']['nz'])
        te = np.array(instate['instate']['te'])
        ti = np.array(instate['instate']['ti'])
        wbeam = np.array(instate['instate']['wbeam'])
        vol = np.array(instate['inmetric']['vol'])

        wth  = 1.5 * 1.602e3 * (ne * te + (ni + nz) * ti)
        betan_th = formula.betan(w=wth, vol=vol, ip=ip, b0=bt, a0=a)
        betan_beam = formula.betan(w=wbeam*1.e6, vol=vol, ip=ip, b0=bt, a0=a)
        betan = betan_th + betan_beam

    elif ps_backend == 'PS':
        ps = plasmastate('ips', 1)
        ps.read(fn_plasma_state)

        rho = ps['rho'][:]
        ne = ps.cell2node(ps['ns'][0, :] * 1.e-19)
        ni = ps.cell2node(np.sum(ps['ns'][1:, :], axis=0) * 1.e-19)
        te = ps.cell2node(ps['Ts'][0, :])
        ti = ps.cell2node(ps['Ti'][:])

        density_beam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0) * 1.e-19
        wbeam = ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0) \
            + ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0)
        wbeam = density_beam*wbeam

        density_alpha = ps.dump_profile(rho, 'rho_fus', 'nfusi', k=0) * 1.e-19
        walpha = ps.dump_profile(rho, 'rho_fus', 'eperp_fusi', k=0) \
            + ps.dump_profile(rho, 'rho_fus', 'epll_fusi', k=0)
        walpha = density_alpha*walpha

        vol = zinterp(ps['rho_eq'], ps['vol'])(rho)
        wtot = 1.5 * 1.602e3 * (ne * te + ni * ti) + 1.602e3*(wbeam + walpha)
        betan = formula.betan(w=wtot, vol=vol, ip=ip, b0=bt, a0=a)

    input_eped = Namelist()
    input_eped['input_eped']['r'] = [r]
    input_eped['input_eped']['a'] = [a]
    input_eped['input_eped']['kappa'] = [kappa]
    input_eped['input_eped']['delta'] = [delta]
    input_eped['input_eped']['bt'] = [bt]
    input_eped['input_eped']['ip'] = [ip]
    input_eped['input_eped']['neped'] = [neped]
    input_eped['input_eped']['nesep'] = [nesep]
    input_eped['input_eped']['zeffped'] = [zeff_ped]
    input_eped['input_eped']['betan'] = [betan]
    input_eped.write('input_eped')


def update_state(fn_inped, fn_instate, fn_plasma_state, ps_backend='INSTATE'):
    wped_model = loglinear(fn_inped, 'wped')
    pped_model = loglinear(fn_inped, 'pped')

    input_eped = Namelist('input_eped')

    wped = wped_model.get_(input_eped['input_eped'])
    pped = pped_model.get_(input_eped['input_eped'])

    ne_ped = input_eped['input_eped']['neped'][0]
    zeff_ped = input_eped['input_eped']['zeffped'][0]

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

    # wped_in = instate['instate']['xwid'][0]
    # rho_ped_in = 1. - wped_in
    # rho_top_in = 1. - 1.5 * wped_in

    if ps_backend == 'INSTATE':
        rho = instate['instate']['rho']
        nrho = instate['instate']['nrho'][0]
        ne = instate['instate']['ne']
        te = instate['instate']['te']
        ti = instate['instate']['ti']
    elif ps_backend == 'PS':
        ps = plasmastate('ips', 1)
        ps.read(fn_plasma_state)
        rho   = ps['rho'][:]
        nrho = len(rho)
        ne = ps.cell2node_bdry(ps['ns'][0, :]) * 1.e-19
        te = ps.cell2node_bdry(ps['Ts'][0, :])
        ti = ps.cell2node_bdry(ps['Ti'][:])
    ne_axis = ne[0]
    te_axis = te[0]
    ti_axis = ti[0]

    ne_sep = instate['instate']['ne_sep'][0]
    te_sep = instate['instate']['te_sep'][0]
    ti_sep = instate['instate']['ti_sep'][0]

    #-- patch updated eped profile
    ne_profile = profile_core_pedestal(x=rho, y=ne, wped=wped, yped=ne_ped, yaxis=ne_axis, ysep=ne_sep)
    te_profile = profile_core_pedestal(x=rho, y=te, wped=wped, yped=te_ped, yaxis=te_axis, ysep=te_sep)
    ti_profile = profile_core_pedestal(x=rho, y=ti, wped=wped, yped=ti_ped, yaxis=ti_axis, ysep=ti_sep)

    ne_update = ne_profile.patch_eped(rho)
    te_update = te_profile.patch_eped(rho)
    ti_update = ti_profile.patch_eped(rho)

    instate['instate']['xwid'] = [wped]
    instate['instate']['xmid'] = [1. - 0.5 * wped]

    # instate['instate']['ne'] = ne_update
    instate['instate']['ne_ped'] = [ne_ped]
    instate['instate']['ne_xwid'] = [wped]
    instate['instate']['ne_xmid'] = [1. - 0.5 * wped]

    # instate['instate']['te'] = te_update
    instate['instate']['te_ped'] = [te_ped]
    instate['instate']['te_xwid'] = [wped]
    instate['instate']['te_xmid'] = [1. - 0.5 * wped]

    # instate['instate']['ti'] = ti_update
    instate['instate']['ti_ped'] = [ti_ped]
    instate['instate']['ti_xwid'] = [wped]
    instate['instate']['ti_xmid'] = [1. - 0.5 * wped]

    if ps_backend == 'INSTATE':
        instate['instate']['ne'] = ne_update
        instate['instate']['te'] = te_update
        instate['instate']['ti'] = ti_update
    elif ps_backend == 'PS':
        ps['ns'][0] = 1.e19 * ps.node2cell(ne_update)
        ps['Ts'][0,:] = ps.node2cell(te_update)
        nspec_th = len(ps['Ts']) - 1
        for k in range(nspec_th):
            ps['Ts'][k + 1, :] = ps.node2cell(ti_update)
        ps['Ti'][:] = ps.node2cell(ti_update)
        ps.store(fn_plasma_state)

    instate.write(fn_instate)
