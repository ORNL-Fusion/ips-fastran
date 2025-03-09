import os
import numpy as np
from Namelist import Namelist
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.timetrace import Timetrace
from fastran.util.model_profile import profile_pedestal
from fastran.util import shape_io
from fastran.util.formula import mu0
from fastran.util.zinterp import zinterp
from scipy.interpolate import interp1d

instate_variables = [
    'ne',
    'te',
    'ti',
    'omega',
    'zeff',
    'j_oh',
    'j_bs',
    'j_nb',
    'j_ec',
    'j_ic',
    'j_lh',
    'pe_nb',
    'pe_ec',
    'pe_ic',
    'pe_lh',
    'pe_fus',
    'pe_ionization',
    'p_rad',
    'pi_nb',
    'pi_ic',
    'pi_lh',
    'pi_fus',
    'pi_ionization',
    'pi_cx',
    'p_ohm',
    'p_ei',
    'torque_nb',
    'torque_in',
    'se_nb',
    'se_ionization',
    'si_nb',
    'si_ionization',
    'q',
    'psipol',
    'psi',
    'density_beam',
    'wbeam',
    'tbeam',
    'density_alpha',
    'walpha',
    'talpha',
    'chie',
    'chii',
    'p_eq',
    'j_tot'
]


class Instate():
    def __init__(self, fname):
        os.system('sync -f ' + fname)
        self.data = Namelist(fname)

        # -- check error
        self.check()

        # -- convert to numpy array
        for key in self.keys():
            if type(self[key][0]) is not type(''):
                self[key] = np.array(self[key])

    def check(self):
        self.density_model = self['density_model'][0]
        if self.density_model not in [0, 1]:
            raise Exception('density_model error')

        self.model_shape = self['model_shape'][0]
        if self.model_shape not in [0, 1, 2]:
            raise Exception('model_shape error')

        self.model_beam = self.default('model_beam', [0])[0]
        if self.model_beam not in [0, 1]:
            raise Exception('model_beam error')
        print(f'model_beam = {self.model_beam}')

    def __getitem__(self, key):
        return self.data['instate'][key]

    def __setitem__(self, key, value):
        self.data['instate'][key] = value

    def __contains__(self, key):
        return True if key.upper() in self.keys() else False

    def keys(self):
        return self.data['instate'].keys()

    def default(self, key, value):
        if key in self:
            return self[key]
        else:
            return value

    def write(self, fname):
        self.data.write(fname)

    def write_netcdf(self, fname):
        pass

    def zeros(self):
        for key in instate_variables:
            if key.upper() not in self.data['instate'].keys():
                self[key] = np.zeros(self['nrho'][0])

    def scale(self):
        for key in self.keys():
            scale = self.default(f'scale_{key}', [-1.])[0]
            if scale > 0.:
                print(f'{key} profile scaled {scale}')
                self[key] = self[key] * scale

    def resize(self, nrho, kind='linear'):
        rho = np.linspace(0, 1, nrho)
        for key in instate_variables:
            if key.upper() in self.data['instate'].keys():
                print(key, len(self['rho']), len(self[key]))
                self[key] = interp1d(self['rho'], self[key], kind=kind, fill_value=(
                    self[key][0], self[key][-1]), bounds_error=False)(rho)
        self['nrho'] = [nrho]
        self['rho'] = rho

    def set_grid(self):
        nrho = self['nrho'][0]
        self['rho'] = np.arange(nrho) / (nrho - 1.)

    def set_shape(self):
        if self.model_shape == 0:
            print('Shape from instate data array')
            rb, zb, rlim, zlim = \
                self['rbdry'], \
                self['zbdry'], \
                self['rlim'], \
                self['zlim']
        elif self.model_shape == 1:
            print('Simple shape')
            rb, zb = \
                shape_io.set_shape(
                    R0=self['r0'][0],
                    a0=self['a0'][0],
                    kappa=self['kappa'][0],
                    delta=self['delta'][0],
                    nt=self['nbdry'][0])
            rlim, zlim = shape_io.set_limiter(rb, zb)
        elif self.model_shape == 2:
            print('Luce shape')
            nbdry = self['nbdry'][0]
            R0 = self['r0'][0]
            a0 = self['a0'][0]
            eps = a0 / R0
            kapu = self['kappa'][0]
            kapl = self['kappa'][0]
            delu = self['delta'][0]
            dell = self['delta'][0]
            z0 = self.default('z0', [0.])[0]
            zetaou = 0.
            zetaiu = 0.
            zetail = 0.
            zetaol = 0.
            zoffset = 0.

            rb, zb, zref = shape_io.boundaryShape(
                a0, eps, kapu, kapl, delu, dell,
                zetaou, zetaiu, zetail, zetaol, zoffset,
                upnull=True, lonull=True, npts=nbdry, doPlot=False)

            if np.isnan(zb[0]):
                zb[0] = 0.0  # need to debug Luce shape routine

            rb = np.append(rb, rb[0])
            zb = np.append(zb, zb[0])

            zb += z0

            rlim, zlim = shape_io.set_limiter(rb, zb)

        self['nbdry'] = [len(rb)]
        self['rbdry'] = rb
        self['zbdry'] = zb
        self['nlim'] = [len(rlim)]
        self['rlim'] = rlim
        self['zlim'] = zlim

    def from_timetrace(
        self,
        f_timetrace,
        timenow,
        interpolation_method='linear',
        init=False,
        force=[]
    ):
        '''
        extract profile timetrace data
        '''
        timetrace = Timetrace(f_timetrace, kind=interpolation_method)

        ip = timetrace.get('ip', timenow) * 1.e-6
        b0 = timetrace.get('bt', timenow)
        r0 = timetrace.get('r0', timenow)
        b0 = np.abs(b0)
        print(f'ip = {ip}, bt={b0}, r0={r0}')

        self['ip'] = [ip]
        self['b0'] = [b0]
        self['r0'] = [r0]

        print(f'time slice at t = {timenow}')
        print('instate from timetrace: update profile')

        for key in ['ne', 'te', 'ti', 'omega', 'nz', 'p_eq', 'j_tot']:
            if self.default(f'trace_{key}', [1])[0] or key in force or init:
                print(f'instate from time trace: {key}')
                self[key] = timetrace.slice(key, timenow, self['rho'])
        if self.default(f'trace_nz', [1])[0] or init:
            print('instate from time trace: nz')
            self['zeff'] = 30. * np.array(self['nz']) / np.array(self['ne']) + 1.
            # self['zeff'][-5:] = self['zeff'][-5] # < ===

        nbdry = timetrace.nearest('nbdry', timenow)
        self['nbdry'] = [nbdry]
        self['rbdry'] = timetrace.nearest('rbdry', timenow)[:nbdry]
        self['zbdry'] = timetrace.nearest('zbdry', timenow)[:nbdry]

        for key in ['p_rad', 'pe_nb', 'pi_nb', 'wbeam', 'density_beam', 'p_ohm']:
            if self.default(f'trace_{key}', [0])[0]:
                print(f'instate from time trace: {key}')
                self[key] = timetrace.slice(key, timenow, self['rho'])

        if self.default(f'trace_nb', [1])[0]:
            for key in ['pnbi']:
                print(f'instate from time trace: {key}')
                self[key] = timetrace.slice_list(key, timenow)

        if self.default(f'trace_ec', [1])[0]:
            for key in ['pech']:
                print(f'instate from time trace: {key}')
                self[key] = timetrace.slice_list(key, timenow)

    def particle_balance(self):
        nrho = self['nrho'][0]
        rho = self['rho']

        n_ion = self['n_ion'][0]
        z_ion = self['z_ion']
        a_ion = self['a_ion']
        f_ion = self['f_ion']

        n_imp = self['n_imp'][0]
        z_imp = self['z_imp']
        a_imp = self['a_imp']
        f_imp = self['f_imp']

        z_beam = self['z_beam'][0]

        density_ion = np.zeros((n_ion, nrho))
        density_imp = np.zeros((n_imp, nrho))
        density_th = np.zeros(nrho)

        if self.density_model == 0:
            print('density_model = 0: impurity density from zeff')

            a = 0
            b = 0
            c = 0
            d = 0
            for k in range(n_imp):
                b = b + f_imp[k] * z_imp[k]
                d = d + f_imp[k] * z_imp[k] * z_imp[k]
            for k in range(n_ion):
                a = a + f_ion[k] * z_ion[k]
                c = c + f_ion[k] * z_ion[k] * z_ion[k]

            for i in range(nrho):
                zne_adj = self['ne'][i]
                zzne_adj = self['ne'][i] * self['zeff'][i]

                # depletion due to beam ions
                zne_adj = zne_adj - z_beam * self['density_beam'][i] - 2. * self['density_alpha'][i]
                zzne_adj = zzne_adj - z_beam**2 * self['density_beam'][i] - 2.**2 * self['density_alpha'][i]

                # effective main ion and impurity densities
                nion = (zne_adj * d - zzne_adj * b) / (a * d - b * c)
                nimp = (zzne_adj * a - zne_adj * c) / (a * d - b * c)

                for k in range(n_ion):
                    density_ion[k][i] = f_ion[k] * nion
                for k in range(n_imp):
                    density_imp[k][i] = f_imp[k] * nimp

        elif self.density_model == 1:
            print('density_model = 1: impurity = f*ne')

            for k in range(n_imp):
                density_imp[k] = self['ne'] * f_imp[k]
                if self['z_imp'][k] == 2 and self['a_imp'][k] == 4 and 'density_he' in self:
                    print('He profile found')
                    density_imp[k] = self['density_he']

            nith = self['ne'] - self['density_beam'] - 2. * self['density_alpha']
            for k in range(n_imp):
                nith = nith - z_imp[k] * density_imp[k]
            for k in range(n_ion):
                density_ion[k] = nith * f_ion[k]

            zeff = nith + z_beam**2 * self['density_beam'] + 2.**2 * self['density_alpha']
            for k in range(n_imp):
                zeff = zeff + z_imp[k]**2 * density_imp[k]
            self['zeff'] = zeff / self['ne']

        for k in range(n_ion):
            self['density_ion_{}'.format(k)] = density_ion[k]
            density_th = density_th + density_ion[k]
        for k in range(n_imp):
            self['density_imp_{}'.format(k)] = density_imp[k]
            density_th = density_th + density_imp[k]
        self['density_th'] = density_th

        self['ni'] = np.array([sum(tmp) for tmp in density_ion.transpose()])
        self['nz'] = np.array([sum(tmp) for tmp in density_imp.transpose()])

    def model_profile(self):
        # -- generate profies from the model profile parameters
        nrho = self['nrho'][0]
        rho = np.arange(nrho) / (nrho - 1.)
        self['rho'] = rho

        # -- default for backward compatibility
        model_profile_list = []
        for key in ['ne', 'te', 'ti']:
            self[f'model_profile_{key}'] = self.default(f'model_profile_{key}', ['ped'])
            model_profile_list.append(key)
        for key in ['density_beam', 'j_tot', 'omega']:
            self[f'model_profile_{key}'] = self.default(f'model_profile_{key}', ['parab'])
            model_profile_list.append(key)
        for key in ['zeff', 'tbeami']:
            self[f'model_profile_{key}'] = self.default(f'model_profile_{key}', ['const'])
            model_profile_list.append(key)
        for key in ['density_alpha', 'walpha']:
            self[f'model_profile_{key}'] = self.default(f'model_profile_{key}', ['zero'])
            model_profile_list.append(key)
        for key in ['wbeam']:
            self[f'model_profile_{key}'] = self.default(f'model_profile_{key}', ['derived'])
            model_profile_list.append(key)
        self['j_tot_axis'] = self.default('jpar_axis', [1.])
        self['j_tot_sep'] = self.default('jpar_sep', [0.])
        self['j_tot_alpha'] = self.default('jpar_alpha', [1.5])
        self['j_tot_beta'] = self.default('jpar_beta', [1.5])
        self['density_beam_axis'] = self.default('nbeam_axis', [0.])
        self['density_beam_sep'] = self.default('nbeam_sep', [0.])
        self['density_beam_alpha'] = self.default('nbeam_alpha', [1.5])
        self['density_beam_beta'] = self.default('nbeam_beta', [1.5])
        self['tbeami_axis'] = self.default('tbeami', [0.])

        xmid = self.default('xmid', [-1.])[0]
        xwid = self.default('xwid', [-1.])[0]
        if xmid > 0:
            self['ne_xmid'] = [xmid]
            self['te_xmid'] = [xmid]
            self['ti_xmid'] = [xmid]
        if xwid > 0:
            self['ne_xwid'] = [xwid]
            self['te_xwid'] = [xwid]
            self['ti_xwid'] = [xwid]

        betan_ped = self.default('betan_ped', [-1.0])[0]
        if betan_ped > 0.0:
            rb = np.array(self['rbdry'])
            zb = np.array(self['zbdry'])
            a0 = 0.5 * (np.max(rb) - np.min(rb))

            b0 = abs(self['b0'][0])
            c_betan = 4.0 * 1.602e5 * self['ne_ped'][0] * mu0 / b0**2 * (a0 * b0 / self['ip'][0])
            te_ped = betan_ped / c_betan
            ti_ped = te_ped
            print('PEDEDSTAL BETAN = ', betan_ped, te_ped)

            self['te_ped'] = [te_ped]
            self['te_xmid'] = [xmid]
            self['te_xwid'] = [xwid]

            self['ti_ped'] = [ti_ped]
            self['ti_xmid'] = [xmid]
            self['ti_xwid'] = [xwid]

        for key in model_profile_list:
            if self[f'model_profile_{key}'][0] == 'ped':
                print(f'{key} from model profile, ped')
                xmid = self[f'{key}_xmid'][0]
                xwid = self[f'{key}_xwid'][0]
                yaxis = self[f'{key}_axis'][0]
                yped = self[f'{key}_ped'][0]
                ysep = self[f'{key}_sep'][0]
                alpha = self[f'{key}_alpha'][0]
                beta = self[f'{key}_beta'][0]
                fit = self[f'{key}_fit'][0]
                self[key] = profile_pedestal(
                    nrho,
                    xmid,
                    xwid,
                    yped,
                    ysep,
                    yaxis,
                    alpha,
                    beta,
                    ifit=fit)(rho)

            elif self[f'model_profile_{key}'][0] == 'parab':
                print(f'{key} from model profile, parab')
                yaxis = self[f'{key}_axis'][0]
                ysep = self[f'{key}_sep'][0]
                alpha = self[f'{key}_alpha'][0]
                beta = self[f'{key}_beta'][0]
                self[key] = (yaxis - ysep) * (1. - rho**alpha)**beta + ysep

            elif self[f'model_profile_{key}'][0] == 'const':
                print(f'{key} from model profile, const')
                yaxis = self[f'{key}_axis'][0]
                self[key] = yaxis * np.ones(nrho)

            elif self[f'model_profile_{key}'][0] == 'zero':
                print(f'{key} from model profile, zero')
                self[key] = np.zeros(nrho)

        if self['model_profile_wbeam'][0] == 'derived':
            print('wbeam from the model profile')
            self['wbeam'] = np.zeros(nrho)
            for i in range(nrho):
                self['wbeam'][i] = 1.5 * 1.602e3 * self['density_beam'][i] * self['tbeami'][i] * 1.e-6

    def to_ps_profile(self, ps, init=False):
        # ne[0] = ne[1]
        # te[0] = te[1]
        # ti[0] = ti[1]
        # omega[0] = omega[1]
        # zeff[0] = zeff[1]

        # -- density
        if self.default('update_ps_ne', [1])[0]:
            ps['ns'][0, :] = 1.e19 * ps.node2cell(self['ne'])
            for k in range(self['n_ion'][0]):
                ps['ns'][k + 1, :] = 1.e19 * ps.node2cell(self['density_ion_{}'.format(k)])
            for k in range(self['n_imp'][0]):
                ps['ns'][k + self['n_ion'][0] + 1, :] = 1.e19 * ps.node2cell(self['density_imp_{}'.format(k)])
            ps['ni'][:] = 1.e19 * ps.node2cell(self['density_th'])

        # -- temperature
        if self.default('update_ps_te', [1])[0]:
            ps['Ts'][0, :] = ps.node2cell(self['te'])

        if self.default('update_ps_ti', [1])[0]:
            for k in range(self['n_ion'][0]):
                ps['Ts'][k + 1, :] = ps.node2cell(self['ti'])
            for k in range(self['n_imp'][0]):
                ps['Ts'][k + self['n_ion'][0] + 1, :] = ps.node2cell(self['ti'])
            ps['Ti'][:] = ps.node2cell(self['ti'])

        # -- zeff
        if self.default('update_nz', [1])[0]:
            ps['Zeff'][:] = ps.node2cell(self['zeff'])
            ps['Zeff_th'][:] = ps.node2cell(self['zeff'])

        # -- rotation
        if self.default('update_omega', [1])[0]:
            ps['omegat'][:] = ps.node2cell(self['omega'])

        # -- radiation
        if self.default('update_p_rad', [1])[0]:
            ps.load_vol_profile(self['rho'], -1.e6 * self['p_rad'], 'rho_rad', 'prad')  # instate p_rad < 0

        # -- ohmic heating (to be removed)
        if self.default('update_p_ohm', [1])[0]:
            ps.load_vol_profile(self['rho'], 1.e6 * self['p_ohm'], 'rho', 'pohme')

        # -- MHD equilibrium current
        if self.default('update_j_tot', [1])[0] or init:
            r0 = self['r0'][0]
            b0 = self['b0'][0]
            self['j_tot'][0] = self['j_tot'][1]
            ps.load_j_parallel(self['rho'], 1.e6 * self['j_tot'], 'rho_eq', 'curt', r0, b0, tot=True)

        # -- MHD equilibrium pressure
        if self.default('update_p_eq', [1])[0] or init:
            ps['P_eq'][:] = self['p_eq']

    def to_ps_nb(self, ps):
        nrho = self['nrho'][0]

        if self.model_beam == 0:
            ps['nbeami'][0][:] = 1.e19 * ps.node2cell(self['density_beam'])
            ps['eperp_beami'][0][:] = 2. * self['tbeam'][0] * np.ones(nrho - 1)
            ps['epll_beami'][0][:] = self['tbeam'][0] * np.ones(nrho - 1)
        elif self.model_beam == 1:
            tbeam = self['wbeam'] / 1.602e-3 / (self['density_beam'] + 1.e-6)
            ps['nbeami'][0][:] = 1.e19 * ps.node2cell(self['density_beam'])
            ps['eperp_beami'][0][:] = 2. / 3. * ps.node2cell(tbeam)
            ps['epll_beami'][0][:] = 1. / 3. * ps.node2cell(tbeam)

        ps.load_vol_profile(self['rho'], 1.e6 * self['pe_nb'], 'rho_nbi', 'pbe')
        ps.load_vol_profile(self['rho'], 1.e6 * self['pi_nb'], 'rho_nbi', 'pbi')
        ps.load_vol_profile(self['rho'], 1.e19 * self['se_nb'], 'rho_nbi', 'sbedep')
        ps.load_vol_profile(self['rho'], self['torque_nb'], 'rho_nbi', 'tqbe')

        r0 = self['r0'][0]
        b0 = self['b0'][0]
        ps.load_j_parallel(self['rho'], 1.e6 * self['j_nb'], 'rho_nbi', 'curbeam', r0, b0)

    def to_ps_rf(self, ps):
        ps.load_vol_profile(self['rho'], 1.e6 * self['pe_ec'], 'rho_ecrf', 'peech')
        ps.load_vol_profile(self['rho'], 1.e6 * self['pe_ic'], 'rho_icrf', 'picrf_totals', k=0)
        ps.load_vol_profile(self['rho'], 1.e6 * self['pi_ic'], 'rho_icrf', 'picrf_totals', k=1)

        r0 = self['r0'][0]
        b0 = self['b0'][0]
        ps.load_j_parallel(self['rho'], 1.e6 * self['j_ec'], 'rho_ecrf', 'curech', r0, b0)
        ps.load_j_parallel(self['rho'], 1.e6 * self['j_ic'], 'rho_icrf', 'curich', r0, b0)

    def to_ps_pnbi(self, ps):
        ps['power_nbi'][:] = self['pnbi']

    def to_ps(self, ps, init=False):
        self.to_ps_profile(ps, init)
        self.to_ps_nb(ps)
        self.to_ps_rf(ps)

    def from_ps(self, ps):
        r0 = self['r0'][0]
        b0 = abs(self['b0'])
        ip = self['ip'][0] * 1.e6

        nrho = len(ps['rho'])
        rho = ps['rho'][:]
        if ps['psipol'][-1] == 0: 
            psi = np.zeros(nrho)
        else:
            psi = ps['psipol'][:] / ps['psipol'][-1]

        ne = ps['ns'][0, :] * 1.e-19
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

        p_eq = ps['P_eq'][:]

        j_tot = ps.dump_j_parallel(rho, 'rho_eq', 'curt', r0, b0, tot=True) * 1.e-6
        j_nb = ps.dump_j_parallel(rho, 'rho_nbi', 'curbeam', r0, b0) * 1.e-6
        j_ec = ps.dump_j_parallel(rho, 'rho_ecrf', 'curech', r0, b0) * 1.e-6
        j_ic = ps.dump_j_parallel(rho, 'rho_icrf', 'curich', r0, b0) * 1.e-6
        j_lh = ps.dump_j_parallel(rho, 'rho_lhrf', 'curlh', r0, b0) * 1.e-6
        j_bs = np.zeros(nrho)
        j_oh = np.zeros(nrho)

        density_beam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0) * 1.e-19
        wbeam = ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0) \
            + ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0)
        wbeam = density_beam * wbeam * 1.602e-3  # MJ/m**3

        pe_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbe') * 1.e-6
        pi_nb = (ps.dump_vol_profile(rho, 'rho_nbi', 'pbi') + ps.dump_vol_profile(rho, 'rho_nbi', 'pbth')) * 1.e-6
        pth_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbth') * 1.e-6

        density_alpha = ps.dump_profile(rho, 'rho_fus', 'nfusi', k=0) * 1.e-19
        walpha = ps.dump_profile(rho, 'rho_fus', 'eperp_fusi', k=0) \
            + ps.dump_profile(rho, 'rho_fus', 'epll_fusi', k=0)
        walpha = density_alpha * walpha * 1.602e-3  # MJ/m**3

        pe_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfuse') * 1.e-6
        pi_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfusi') * 1.e-6
        pth_fus = ps.dump_vol_profile(rho, 'rho_fus', 'pfusth') * 1.e-6

        pe_ec = ps.dump_vol_profile(rho, 'rho_ecrf', 'peech') * 1.e-6
        pe_ic = ps.dump_vol_profile(rho, 'rho_icrf', 'picrf_totals', k=0) * 1.e-6
        pe_lh = ps.dump_vol_profile(rho, 'rho_lhrf', 'pelh') * 1.e-6

        pi_ic = ps.dump_vol_profile(rho, 'rho_icrf', 'picrf_totals', k=1) * 1.e-6
        pi_lh = ps.dump_vol_profile(rho, 'rho_lhrf', 'pilh') * 1.e-6

        tqbe = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbe')
        tqbi = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbi')
        tqbjxb = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbjxb')
        tqbth = ps.dump_vol_profile(rho, 'rho_nbi', 'tqbth')

        torque_nb = tqbe + tqbi + tqbjxb + tqbth
        torque_in = np.zeros(nrho)

        se_nb = (ps.dump_vol_profile(rho, 'rho_nbi', 'sbedep') + ps.dump_vol_profile(rho, 'rho_nbi', 'sbehalo')) * 1.e-19

        p_ei = np.zeros(nrho)
        p_rad = ps.dump_vol_profile(rho, 'rho_rad', 'prad') * 1.e-6
        p_ohm = ps.dump_vol_profile(rho, 'rho_rad', 'pohme') * 1.e-6
        pe_ionization = np.zeros(nrho)
        pi_ionization = np.zeros(nrho)
        pi_cx = np.zeros(nrho)
        si_nb = np.zeros(nrho)
        chie = np.zeros(nrho)
        chii = np.zeros(nrho)

        se_ionization = np.zeros(nrho)
        si_ionization = np.zeros(nrho)

        pe_trans = ps.dump_vol_profile(rho, 'rho', 'pe_trans') * 1.e-6
        pi_trans = ps.dump_vol_profile(rho, 'rho', 'pi_trans') * 1.e-6

        self['nrho'] = [nrho]
        self['rho'] = rho
        self['ne'] = ne
        self['te'] = te
        self['ti'] = ti
        self['zeff'] = zeff
        self['omega'] = omega
        self['p_eq'] = p_eq
        self['j_tot'] = j_tot
        self['j_oh'] = j_oh
        self['j_bs'] = j_bs
        self['j_nb'] = j_nb
        self['j_ec'] = j_ec
        self['j_ic'] = j_ic
        self['j_lh'] = j_lh
        self['pe_nb'] = pe_nb
        self['pe_ec'] = pe_ec
        self['pe_ic'] = pe_ic
        self['pe_lh'] = pe_lh
        self['pe_fus'] = pe_fus
        self['pe_ionization'] = pe_ionization
        self['p_rad'] = -p_rad  # instate p_rad < 0
        self['pi_nb'] = pi_nb
        self['pi_ic'] = pi_ic
        self['pi_lh'] = pi_lh
        self['pi_fus'] = pi_fus
        self['pi_ionization'] = pi_ionization
        self['pi_cx'] = pi_cx
        self['p_ohm'] = p_ohm
        self['p_ei'] = p_ei
        self['torque_nb'] = torque_nb
        self['torque_in'] = torque_in
        self['se_nb'] = se_nb
        self['se_ionization'] = se_ionization
        self['si_nb'] = si_nb
        self['si_ionization'] = si_ionization
        self['q'] = q
        self['psi'] = psi
        self['density_beam'] = density_beam
        self['wbeam'] = wbeam
        self['density_alpha'] = density_alpha
        self['walpha'] = walpha
        self['chie'] = chie
        self['chii'] = chii

        self['pe_trans'] = pe_trans
        self['pi_trans'] = pi_trans

    def inmetric(self, f_geqdsk, nrho=101, nth=101):
        pi = np.pi

        r0 = self['r0'][0]
        b0 = abs(self['b0'][0])
        ip = self['ip'][0]
        ps = plasmastate('ips', 1)
        ps.init_from_geqdsk(f_geqdsk, nrho=int(nrho), nth=int(nth))

        psi = ps['psipol'][:] / ps['psipol'][-1]
        rho = np.sqrt(ps['phit'][:] / ps['phit'][-1])
        nrho = len(rho)
        rhob = (ps['phit'][-1] / pi / abs(b0))**0.5
        rhopsi = zinterp(psi, rho)
        ipol = ps['g_eq'][:] / (r0 * b0)
        ipol = abs(ipol)
        volp = 4.0 * pi * pi * rho * rhob * r0 / ipol / (r0 * r0 * ps['gr2i'][:])
        g11 = volp * ps['grho2'][:] * rhob**2
        g22 = r0 * volp / (4.0 * pi * pi) * ps['grho2r2i'][:] * rhob**2
        g33 = r0 * r0 * ps['gr2i'][:]
        gradrho = ps['grho1'][:] * rhob
        area = ps['surf'][:]
        rmajor = ps['Rmajor_mean'][:]
        rminor = ps['rMinor_mean'][:]
        shift = rmajor - r0
        kappa = ps['elong'][:]
        delta = ps['triang'][:]
        pmhd = ps['P_eq'][:]
        qmhd = ps['q_eq'][:]
        nc1 = ps['gncfh'][:]
        gb1 = ps['gb1'][:]
        gb2 = ps['gb2'][:]
        Bmax = ps['B_surfMax'][:]
        hfac1 = gb1 / Bmax
        hfac2 = gb2 / Bmax**2
        bp2 = ps['gb2'][:] - ps['g_eq'][:]**2 * ps['gr2i'][:]

        self.data['inmetric']['ip'] = [ip * 1.e6]
        self.data['inmetric']['bcentr'] = [b0]
        self.data['inmetric']['rmajor'] = [r0]
        self.data['inmetric']['aminor'] = [rminor[-1]]
        self.data['inmetric']['kappa'] = [kappa[-1]]
        self.data['inmetric']['delta'] = [delta[-1]]
        self.data['inmetric']['nrho'] = [nrho]
        self.data['inmetric']['rhob'] = [rhob]
        self.data['inmetric']['rho'] = rho
        self.data['inmetric']['volp'] = volp
        self.data['inmetric']['ipol'] = ipol
        self.data['inmetric']['g11'] = g11
        self.data['inmetric']['g22'] = g22
        self.data['inmetric']['g33'] = g33
        self.data['inmetric']['gradrho'] = gradrho
        self.data['inmetric']['area'] = area
        self.data['inmetric']['rminor'] = rminor
        self.data['inmetric']['rmajor'] = rmajor
        self.data['inmetric']['shift'] = shift
        self.data['inmetric']['kappa'] = kappa
        self.data['inmetric']['delta'] = delta
        self.data['inmetric']['pmhd'] = pmhd
        self.data['inmetric']['qmhd'] = qmhd
        self.data['inmetric']['er'] = np.zeros(nrho)  # <======
        self.data['inmetric']['nc1'] = nc1
        self.data['inmetric']['hfac1'] = hfac1
        self.data['inmetric']['hfac2'] = hfac2
        self.data['inmetric']['psi'] = ps['psipol'][:]
        self.data['inmetric']['vol'] = ps['vol'][:]
        self.data['inmetric']['gr2i'] = ps['gr2i'][:]
        self.data['inmetric']['bp2'] = bp2


if __name__ == '__main__':
    instate = Instate('instate0')

    instate.model_profile()
    instate.zeros()
    instate.particle_balance()
    print(instate.data)

    ps = plasmastate('ips', 1)

    ps.init(
        ps,
        global_label=['ips'],
        runid=['ips'],
        shot_number=[123456],
        nspec_th=[instate['n_ion'][0] + instate['n_imp'][0]],
        nspec_beam=[instate['n_beam'][0]],
        nspec_fusion=[instate['n_fusion'][0]],
        nspec_rfmin=[instate['n_min'][0]],
        nspec_gas=[instate['n_ion'][0]],
        z_ion=instate['z_ion'] + instate['z_imp'],
        a_ion=instate['a_ion'] + instate['a_imp'],
        z_beam=instate['z_beam'],
        a_beam=instate['a_beam'],
        z_fusion=instate['z_fusion'],
        a_fusion=instate['a_fusion'],
        z_rfmin=instate['z_min'],
        a_rfmin=instate['a_min'],
        z_gas=instate['z_ion'],
        a_gas=instate['a_ion'],
        nicrf_src=[1],
        nlhrf_src=[1],
        nrho=[101],
        time=[0.0],
        nlim=instate['nlim']
    )

    instate.to_ps(ps)

    ps.store('ps.nc')
