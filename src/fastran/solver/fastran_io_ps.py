"""
 ----------------------------------------------------------------------
 fastran solver component IO, plasma state backend
 ----------------------------------------------------------------------
"""

import os
import numpy as np
import netCDF4
from Namelist import Namelist
from fastran.util.zinterp import zinterp
from fastran.solver import zfdat
from fastran.equilibrium.efit_eqdsk import readg
from fastran.plasmastate.plasmastate import plasmastate


def write_input(f_state, f_eqdsk, rdir='.', recycle=0.):
    # -- read geqdsk and plasma state file
    geq = readg(f_eqdsk)
    r0 = geq['rzero']
    b0 = np.abs(geq['bcentr'])
    ip = geq['cpasma']

    ps = plasmastate('ips', 1)
    ps.read(f_state)
    spec = ps.get_species()

    # --  profiles
    nrho = len(ps['rho'])
    rho = ps['rho'][:]
    ne = ps['ns'][0, :]*1.e-19
    te = ps['Ts'][0, :]
    ti = ps['Ti'][:]
    zeff = ps['Zeff'][:]
    omega = ps['omegat'][:]
    ne = ps.cell2node_bdry(ne)
    te = ps.cell2node_bdry(te)
    ti = ps.cell2node_bdry(ti)
    zeff = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)
    q = np.zeros(nrho)
    fp = np.zeros(nrho)

    j_tot = ps.dump_j_parallel(rho, 'rho_eq', 'curt', r0, b0, tot=True)*1.e-6
    j_nb = ps.dump_j_parallel(rho, 'rho_nbi', 'curbeam', r0, b0)*1.e-6
    j_rf = ps.dump_j_parallel(rho, 'rho_ecrf', 'curech', r0, b0) \
        + ps.dump_j_parallel(rho, 'rho_icrf', 'curich', r0, b0) \
        + ps.dump_j_parallel(rho, 'rho_lhrf', 'curlh', r0, b0)
    j_rf *= 1.e-6
    j_bs = np.zeros(nrho)

    density_beam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0)*1.e-19
    wbeam = ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0) \
        + ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0)
    wbeam = density_beam*wbeam*1.602e-3 # MJ/m**3

    pe_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbe')*1.e-6
    pi_nb = (ps.dump_vol_profile(rho, 'rho_nbi', 'pbi') + ps.dump_vol_profile(rho, 'rho_nbi', 'pbth'))*1.e-6
    pth_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbth')*1.e-6

    density_alpha = ps.dump_profile(rho, 'rho_fus', 'nfusi', k=0)*1.e-19
    walpha = ps.dump_profile(rho, 'rho_fus', 'eperp_fusi', k=0) \
        + ps.dump_profile(rho, 'rho_fus', 'epll_fusi', k=0)
    walpha = density_alpha*walpha*1.602e-3

    pe_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfuse')*1.e-6
    pi_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfusi')*1.e-6
    pth_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfusth')*1.e-6

    pe_rf = ps.dump_vol_profile(rho, 'rho_ecrf', 'peech') \
        + ps.dump_vol_profile(rho, 'rho_icrf', 'picrf_totals', k=0) \
        + ps.dump_vol_profile(rho, 'rho_lhrf', 'pelh')
    pe_rf *= 1.e-6

    pi_rf = ps.dump_vol_profile(rho, 'rho_icrf', 'picrf_totals', k=1) \
        + ps.dump_vol_profile(rho, 'rho_lhrf', 'pilh')
    pi_rf *= 1.e-6

    tqbe = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbe')
    tqbi = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbi')
    tqbjxb = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbjxb')
    tqbth = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbth')

    torque_nb = tqbe + tqbi + tqbjxb + tqbth
    torque_in = np.zeros(nrho)

    se_nb = (ps.dump_vol_profile(rho, 'rho_nbi', 'sbedep') + ps.dump_vol_profile(rho, 'rho_nbi', 'sbehalo'))*1.e-19

    se_tot = ps.integrate('sbedep')
    print('se_tot =', se_tot)

    if recycle > 0:
        se_recycle = 1.e-19*ps.dump_vol_profile(rho, 'rho_gas', 'sprof0e', k=0)
        se_recycle_tot = ps.integrate('sprof0e', k=0)
        print('se_recycle_tot =', se_recycle_tot)
        se_recycle *= recycle*se_tot
    else:
        se_recycle = np.zeros(nrho)
        print('se_recycle_tot = 0')

    se_pellet = ps.dump_vol_profile(rho, 'rho', 'sn_trans', k=0)*1.e-19

    # p_rad = np.zeros(nrho)
    # p_ohm = np.zeros(nrho)
    p_rad = ps.dump_vol_profile(rho, 'rho_rad', 'prad')*1.e-6
    p_ohm = ps.dump_vol_profile(rho, 'rho', 'pohme')*1.e-6

    pe_ionization = np.zeros(nrho)
    pi_ionization = np.zeros(nrho)
    pi_cx = np.zeros(nrho)
    si_nb = np.zeros(nrho)
    chie = np.zeros(nrho)
    chii = np.zeros(nrho)

    te[-1] = np.abs(te[-1])
    ti[-1] = np.abs(ti[-1])

    amain = spec['amain']
    zmain = spec['zmain']
    n_imp = spec['n_imp']
    a_imp = spec['a_imp']
    z_imp = spec['z_imp']
    f_imp = spec['f_imp']
    f_imp_0 = [f_imp[k][0] for k in range(n_imp)]

    # -- metrics
    psi = ps['psipol'][:]/ps['psipol'][-1]  # equi-drho grid
    rho = np.sqrt(ps['phit'][:]/ps['phit'][-1])
    rhob = (ps['phit'][-1]/np.pi/b0)**0.5

    rhopsi = zinterp(psi, rho)
    ipol = ps['g_eq'][:]/(r0*b0)
    ipol = np.abs(ipol)
    volp = 4.0*np.pi*np.pi*rho*rhob*r0/ipol/(r0*r0*ps['gr2i'][:])

    g11 = volp*ps['grho2'][:]*rhob**2
    g22 = r0*volp/(4.0*np.pi*np.pi)*ps['grho2r2i'][:]*rhob**2
    g33 = r0*r0*ps['gr2i'][:]
    gradrho = ps['grho1'][:]*rhob
    area = ps['surf'][:]
    rmajor = ps['Rmajor_mean'][:]
    rminor = ps['rMinor_mean'][:]
    shift = rmajor-r0
    kappa = ps['elong'][:]
    delta = ps['triang'][:]
    pmhd = ps['P_eq'][:]
    qmhd = ps['q_eq'][:]

    nc1 = ps['gncfh'][:]
    gb1 = ps['gb1'][:]
    gb2 = ps['gb2'][:]
    Bmax = ps['B_surfMax'][:]
    hfac1 = gb1/Bmax
    hfac2 = gb2/Bmax**2

    # -- write inprof
    f = open(os.path.join(rdir, 'inprof'), 'w')
    f.write(' # generated by zfastran.py Ver+Oct2009\n')
    zfdat.write_f('inflag', [1.0], '', f)
    zfdat.write_f('time0', [0.0], 's', f)
    zfdat.write_f('ip', [geq['cpasma']*1.0e-6], '', f)
    zfdat.write_f('bcentr', [geq["bcentr"]], '', f)
    zfdat.write_f('rmajor', [r0], '', f)
    zfdat.write_f('aminor', [rminor[-1]], '', f)
    zfdat.write_f('elongb', [kappa[-1]], '', f)
    zfdat.write_f('trianb', [delta[-1]], '', f)
    zfdat.write_f('rho', rho, '', f)
    zfdat.write_f('amain', [amain[0]], '', f)
    zfdat.write_f('zmain', [zmain[0]], '', f)
    zfdat.write_i('nimp', [n_imp], '', f)
    zfdat.write_f('aimp', a_imp, '', f)
    zfdat.write_f('zimp', z_imp, '', f)
    zfdat.write_f('fimp', f_imp_0, '', f)
    zfdat.write_f('ene', ne, '', f)
    zfdat.write_f('te', te, '', f)
    zfdat.write_f('ti', ti, '', f)
    zfdat.write_f('zeff', zeff, '', f)
    zfdat.write_f('omega', omega, '', f)
    zfdat.write_f('q', q, '', f)
    zfdat.write_f('fp', fp, '', f)
    zfdat.write_f('curpar', j_tot, '', f)
    zfdat.write_f('curbeam', j_nb, '', f)
    zfdat.write_f('currf', j_rf, '', f)
    zfdat.write_f('curboot', j_bs, '', f)
    zfdat.write_f('enbeam', density_beam, '', f)
    zfdat.write_f('wbeam', wbeam, '', f)
    zfdat.write_f('qbeame', pe_nb, '', f)
    zfdat.write_f('qbeami', pi_nb, '', f)
    zfdat.write_f('enalp', density_alpha, '', f)
    zfdat.write_f('walp', walpha, '', f)
    zfdat.write_f('qtfuse', pe_fus, '', f)
    zfdat.write_f('qtfusi', pi_fus, '', f)
    zfdat.write_f('qrfe', pe_rf, '', f)
    zfdat.write_f('qrfi', pi_rf, '', f)
    zfdat.write_f('qrad', p_rad, '', f)
    zfdat.write_f('qohm', p_ohm, '', f)
    zfdat.write_f('qione', pe_ionization, '', f)
    zfdat.write_f('qioni', pi_ionization, '', f)
    zfdat.write_f('qcx', pi_cx, '', f)
    zfdat.write_f('storqueb', torque_nb, '', f)
    zfdat.write_f('storque', torque_in, '', f)
    zfdat.write_f('sion', se_nb + se_recycle + se_pellet, '', f)
    zfdat.write_f('chie', chie, '', f)
    zfdat.write_f('chii', chii, '', f)
    f.close()

    # -- write inmetric
    f = open(os.path.join(rdir, 'inmetric'), 'w')
    zfdat.write_f('Ip', [ip], '', f)
    zfdat.write_f('bcentr', [b0], '', f)
    zfdat.write_f('rmajor', [r0], '', f)
    zfdat.write_f('aminor', [rminor[-1]], '', f)
    zfdat.write_f('kappa', [kappa[-1]], '', f)
    zfdat.write_f('delta', [delta[-1]], '', f)
    zfdat.write_i('nrho', [nrho], '', f)
    zfdat.write_f('rhob', [rhob], '', f)
    zfdat.write_f('rho', rho, '', f)
    zfdat.write_f('volp', volp, '', f)
    zfdat.write_f('ipol', ipol, '', f)
    zfdat.write_f('g11', g11, '', f)
    zfdat.write_f('g22', g22, '', f)
    zfdat.write_f('g33', g33, '', f)
    zfdat.write_f('gradrho', gradrho, '', f)
    zfdat.write_f('area', area, '', f)
    zfdat.write_f('a', rminor, '', f)
    zfdat.write_f('rtor', rmajor, '', f)
    zfdat.write_f('shift', shift, '', f)
    zfdat.write_f('elong', kappa, '', f)
    zfdat.write_f('triag', delta, '', f)
    zfdat.write_f('pmhd', pmhd, '', f)
    zfdat.write_f('qmhd', qmhd, '', f)
    zfdat.write_f('er', np.zeros(nrho), '', f)
    zfdat.write_f('nc1', nc1, '', f)
    zfdat.write_f('hfac1', hfac1, '', f)
    zfdat.write_f('hfac2', hfac2, '', f)
    f.close()


def update_state(f_state, f_eqdsk, f_fastran, f_instate, relax=1., relax_j=1., adjust_ip=0, fni_target=1., relax_ip=0.3):
    # -- read geqdsk and plasma state
    ps = plasmastate('ips', 1)
    ps.read(f_state)
    spec = ps.get_species()

    geq = readg(f_eqdsk)
    r0 = geq['rzero']
    b0 = np.abs(geq['bcentr'])
    ip = geq['cpasma']

    # -- read fastran
    fastran = netCDF4.Dataset(f_fastran, 'r', format='NETCDF4')
    rho = fastran.variables['rho'][:]
    nrho = len(rho)

    if nrho != len(ps['rho']):
        raise Exception('nrho differ, fastran: %d, ps: %s' % (nrho, ps.nrho))

    ne = fastran.variables['ne'][-1, :]
    te = fastran.variables['te'][-1, :]
    ti = fastran.variables['ti'][-1, :]
    omega = fastran.variables['omega'][-1, :]

    j_tot = fastran.variables['j_tot'][-1, :]*1.e6
    j_oh = fastran.variables['j_oh'][-1, :]*1.e6
    j_bs = fastran.variables['j_bs'][-1, :]*1.e6

    zeff = fastran.variables['zeff'][-1, :]
    ni = fastran.variables['ni'][-1, :]
    nhe = fastran.variables['nhe'][-1, :]
    z_imp = fastran.variables['zimp'][:]
    n_imp = len(z_imp)
    nz = n_imp*[0]
    for k in range(n_imp):
        nz[k] = fastran.variables['nz%d' % k][-1, :]
    nz = np.array(nz)


    # -- update plasma state
    # ps['ns'][0] = (1. - relax)*ps['ns'][0] + 1.e19*relax*ps.node2cell(ne)
    ps['ns'][0] = 1.e19*ps.node2cell(ne)
    for k in range(spec['n_ion']):
        ps['ns'][spec['k_ion'][k], :] = 1.e19*ps.node2cell(ni*spec['f_ion'][k])
    for k in range(spec['n_imp']):
        ps['ns'][spec['k_imp'][k], :] = 1.e19*ps.node2cell(nz[k])
    if spec['k_he4'] > 0:
        ps['ns'][spec['k_he4'], :] = 1.e19*ps.node2cell(nhe)
    ps['Zeff'][:] = ps.node2cell(zeff)
    ps['Zeff_th'][:] = ps.node2cell(zeff)

    ps['Ts'][0, :] = (1. - relax)*ps['Ts'][0, :] + relax*ps.node2cell(te)
    nspec_th = len(ps['Ts']) - 1
    print('nspec_th =', nspec_th)
    for k in range(nspec_th):
        ps['Ts'][k + 1, :] = (1. - relax)*ps['Ts'][k + 1, :] + relax*ps.node2cell(ti)
    ps['Ti'][:] = (1. - relax)*ps['Ti'][:] + relax*ps.node2cell(ti)

    ps['omegat'][:] = (1. - relax)*ps['omegat'][:] + relax*ps.node2cell(omega)

    ps.load_j_parallel(rho, j_tot, 'rho_eq', 'curt', r0, b0, tot=True)
    ps.load_j_parallel(rho, j_bs, 'rho', 'curr_bootstrap', r0, b0)
    ps.load_j_parallel(rho, j_oh, 'rho', 'curr_ohmic', r0, b0)

    pe_fus = fastran.variables['pe_fus'][-1, :]*1.e6
    pi_fus = fastran.variables['pi_fus'][-1, :]*1.e6

    ps.load_vol_profile(rho, pe_fus, 'rho_fus', 'pfuse')
    ps.load_vol_profile(rho, pi_fus, 'rho_fus', 'pfusi')

    ps.update_particle_balance()

    #-- write plasma state
    ps.store(f_state)

    # -- global parameters to instate
    instate = Namelist(f_instate)

    instate['instate']['pnbe'] = [fastran.variables['pnbe'][:][-1]]
    instate['instate']['pnbi'] = [fastran.variables['pnbi'][:][-1]]
    instate['instate']['prfe'] = [fastran.variables['prfe'][:][-1]]
    instate['instate']['prfi'] = [fastran.variables['prfi'][:][-1]]
    instate['instate']['pfuse'] = [fastran.variables['pfuse'][:][-1]]
    instate['instate']['pfusi'] = [fastran.variables['pfusi'][:][-1]]
    instate['instate']['pei'] = [fastran.variables['pei'][:][-1]]
    instate['instate']['prad'] = [fastran.variables['prad'][:][-1]]
    instate['instate']['area'] = [fastran.variables['area'][-1, :][-1]]

    if adjust_ip > 0:
        ibs = fastran.variables['ibs'][-1]
        inb = fastran.variables['inb'][-1]
        irf = fastran.variables['irf'][-1]
        fni = (ibs + inb + irf)/(ip*1.e-6)
        ip1 = (ibs + inb + irf)*fni_target

        ip0 = instate['instate']['ip'][0]
        instate['instate']['ip'][0] = (1. - relax_ip)*ip0 + relax_ip*ip1

        # j_tot = (1. - relax_ip)*j_tot + relax_ip*j_tot/fni

        print('******* IP ADJUST (OLD, NEW):', ip0, ip1)

    # j_tot_prev = ps.dump_j_parallel(rho, 'rho_eq', 'curt', r0, b0, tot=True)
    # j_tot = relax_j*j_tot + (1. - relax_j)*j_tot_prev

    instate.write(f_instate)

# -----------------------------------------------------------------------
# test
if __name__ == "__main__":
    pass
