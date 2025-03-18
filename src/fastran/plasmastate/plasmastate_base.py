"""
 -----------------------------------------------------------------------
 plasma state utility, base api
 -----------------------------------------------------------------------
"""

import numpy as np
from fastran.util.zinterp import zinterp as interp
from Namelist import Namelist


class plasmastate_base():
    def __init__(self):
        pass

    def analytic_volume(self, b0, r0, a0, kappa0, delta0):
        nrho_eq = 101  # self['nrho_eq']
        print('usng analytic volume')

        for i in range(nrho_eq):
            rho = self['rho_eq'][i]
            a = a0 * rho
            kappa = 0.5 * (1. + kappa0 + (kappa0 - 1.) * rho**2)
            delta = delta0 * rho**2

            self['vol'][i] = np.pi * 2. * np.pi * a**2 * \
                kappa * (r0 - 0.25 * a * kappa)
            self['g_eq'][i] = r0 * abs(b0)

        out = Namelist()
        out['check']['vol'] = self['vol'][:]
        out['check']['g_eq'] = self['g_eq'][:]
        out.write('vol.dat')

    def load_vol_profile(self, rho, prof, xkey, ykey, k=-1, sum=0, add=0):
        zone_spl = interp(self['rho_eq'], self['vol'])
        prof_spl = interp(rho, prof)

        rho_ps = self[xkey]
        dat = np.zeros(len(rho_ps) - 1)
        vol_integral = 0.
        for i in range(len(rho_ps) - 1):
            dat[i] = 0.5 * (prof_spl(rho_ps[i + 1]) + prof_spl(rho_ps[i])) \
                * (zone_spl(rho_ps[i + 1]) - zone_spl(rho_ps[i]))
            vol_integral += dat[i]
        if sum > 0:
            dat *= sum / vol_integral

        print(xkey, ykey, len(self[ykey]), len(dat))
        if k < 0:
            if add > 0:
                self[ykey][:] = self[ykey][:] + dat
            else:
                self[ykey][:] = dat
        else:
            if add > 0:
                self[ykey][k, :] = self[ykey][k, :] + dat
            else:
                self[ykey][k, :] = dat

        print('integrated : ', vol_integral)

    def load_vol_profile_m(self, rho, prof, xkey, ykey, ksrc, k):
        zone_spl = interp(self['rho_eq'], self['vol'])
        prof_spl = interp(rho, prof)

        rho_ps = self[xkey]
        dat = np.zeros(len(rho_ps) - 1)
        vol_integral = 0.
        for i in range(len(rho_ps) - 1):
            dat[i] = 0.5 * (prof_spl(rho_ps[i + 1]) + prof_spl(rho_ps[i])) \
                * (zone_spl(rho_ps[i + 1]) - zone_spl(rho_ps[i]))
            vol_integral += dat[i]

        self[ykey][k, ksrc][:] = self[ykey][:]

        print('integrated : ', vol_integral)

    def dump_vol_profile(self, rho, xkey, ykey, k=-1):
        if xkey not in self or ykey not in self:
            print(f'{xkey} or {ykey} not in plasma state')
            return np.zeros(len(rho))

        zone_spl = interp(self['rho_eq'], self['vol'])

        rho_ps = self[xkey]

        if k < 0:
            prof_ps = self[ykey]
        else:
            prof_ps = self[ykey][k]

        cell = np.zeros(len(rho_ps) - 1)
        sum = 0.
        for i in range(len(rho_ps) - 1):
            cell[i] = prof_ps[i] / \
                (zone_spl(rho_ps[i + 1]) - zone_spl(rho_ps[i]))
            sum += prof_ps[i]
        node = self.cell2node(cell)

        print('dump_vol_profile :', xkey, ykey, sum)

        return interp(rho_ps, node)(rho)

    def dump_vol_profile_m(self, rho, xkey, ykey, ksrc, k):
        if xkey not in self or ykey not in self:
            print(f'{xkey} or {ykey} not in plasma state')
            return np.zeros(len(rho))

        zone_spl = interp(self['rho_eq'], self['vol'])

        rho_ps = self[xkey]
        prof_ps = self[ykey][ksrc, k]

        cell = np.zeros(len(rho_ps) - 1)
        sum = 0.
        for i in range(len(rho_ps) - 1):
            cell[i] = prof_ps[i] / \
                (zone_spl(rho_ps[i + 1]) - zone_spl(rho_ps[i]))
            sum += prof_ps[i]
        node = self.cell2node(cell)

        print('dump_vol_profile :', xkey, ykey, sum)

        return interp(rho_ps, node)(rho)

    def integrate(self, ykey, k=-1):
        if k < 0:
            prof_ps = self[ykey]
        else:
            prof_ps = self[ykey][k]

        sum = 0.
        for i in range(len(prof_ps)):
            sum += prof_ps[i]
        return sum

    def load_profile(self, rho, prof, xkey, ykey, k=-1):
        prof_spl = interp(rho, prof)

        rho_ps = self[xkey]
        prof_ps = prof_spl(rho_ps)
        prof_ps = self.node2cell(prof_ps)
        if k < 0:
            self[ykey][:] = prof_ps
        else:
            self[ykey][k, :] = prof_ps

    def dump_profile(self, rho, xkey, ykey, k=-1):
        if xkey not in self or ykey not in self:
            print(f'{xkey} or {ykey} not in plasma state')
            return np.zeros(len(rho))

        rho_ps = self[xkey]

        if k < 0:
            prof_ps = self[ykey]
        else:
            prof_ps = self[ykey][k]
        node = self.cell2node(prof_ps)

        return interp(rho_ps, node)(rho)

    # 0---1----2----3----
    #   0    1    2

    def cell2node(self, cell):
        nrho = len(cell) + 1
        node = np.zeros(nrho)
        node[0] = cell[0]
        for i in range(1, nrho - 1):
            node[i] = 0.5 * (cell[i - 1] + cell[i])
        node[-1] = cell[-1]
        return node

    def node2cell(self, node):
        nrho = len(node)
        cell = np.zeros(nrho - 1)
        cell[0] = node[1]
        for i in range(1, nrho - 1):
            cell[i] = 0.5 * (node[i] + node[i + 1])
        return cell

    def cell2node_bdry(self, cell):
        nrho = len(cell) + 1
        node = np.zeros(nrho)
        node[0] = cell[0]
        node[1] = cell[0]
        for k in range(2, nrho):
            node[k] = 2. * cell[k - 1] - node[k - 1]
        return node

    def load_j_parallel(
            self,
            rho,
            jpar,
            xkey,
            ykey,
            r0,
            b0,
            tot=False,
            add=0,
            ksrc=-1):
        #       if k>=0:
        rho_ps = self[xkey][:]

        ipol = interp(self['rho_eq'], self['g_eq'] / (r0 * b0))(rho_ps)
        vol = interp(self['rho_eq'], self['vol'])(rho_ps)
        area = interp(self['rho_eq'], self['area'])(rho_ps)

        jpar_ps = interp(rho, jpar)(rho_ps)

        curt = np.zeros(len(rho_ps))
        curt[0] = 0.
        for i in range(len(rho_ps) - 1):
            dV = vol[i + 1] - vol[i]
            jparm = 0.5 * (jpar_ps[i] + jpar_ps[i + 1])
            ipolm = 0.5 * (ipol[i] + ipol[i + 1])
            curt[i + 1] = (curt[i] / ipol[i] + jparm * dV /
                           (2. * np.pi * r0 * ipolm**2)) * ipol[i + 1]

        print('I%s=%5.3e' % (ykey, 1.e-6 * curt[-1]) + 'MA')

        if tot:
            self[ykey][:] = curt
        else:
            if add > 0:
                for i in range(len(rho) - 1):
                    self[ykey][i] += curt[i + 1] - curt[i]
            else:
                for i in range(len(rho) - 1):
                    self[ykey][i] = curt[i + 1] - curt[i]

    def dump_j_parallel(self, rho, xkey, ykey, r0, b0, tot=False, k=-1):
        if xkey not in self or ykey not in self:
            print(f'{xkey} or {ykey} not in plasma state')
            return np.zeros(len(rho))

        rho_ps = self[xkey]
        jpar_ps = self[ykey]

        ipol = interp(self['rho_eq'], self['g_eq'][:] / (r0 * b0))(rho_ps)
        vol = interp(self['rho_eq'], self['vol'][:])(rho_ps)
        area = interp(self['rho_eq'], self['area'][:])(rho_ps)

        if tot:
            curt = jpar_ps
        else:
            curt = np.zeros(len(rho_ps))
            curt[0] = 0.
            for i in range(len(rho) - 1):
                curt[i + 1] = curt[i] + jpar_ps[i]

        jpar = np.zeros(len(rho_ps) - 1)
        for i in range(len(rho_ps) - 1):
            dV = vol[i + 1] - vol[i]
            ipolm = 0.5 * (ipol[i] + ipol[i + 1])
            jpar[i] = 2. * np.pi * r0 * ipolm**2 / dV * \
                (curt[i + 1] / ipol[i + 1] - curt[i] / ipol[i])

        return interp(rho, self.cell2node(jpar))(rho)

    def get_species(self):
        rho = self['rho']
        nrho = len(rho)

        ps_xe = 1.6022e-19
        ps_mp = 1.6726e-27

        z_spec = [round(x) for x in self['qatom_S'][1:] / ps_xe]
        a_spec = [round(x) for x in self['m_S'][1:] / ps_mp]
        n_spec = len(z_spec)

        z_ion = []
        a_ion = []
        k_ion = []
        z_imp = []
        a_imp = []
        k_imp = []
        ni = []
        nz = []
        f_imp = []
        nhe4 = np.zeros(nrho)
        k_he4 = -1
        # print self['S_name'][1:]
        print(z_spec)
        spec_list = [''.join(key).strip() for key in self['S_name'][1:]]
        print(spec_list)

        for k in range(n_spec):
            if spec_list[k] in ['H', 'D', 'T']:
                print('ion', spec_list[k])
                z_ion.append(z_spec[k])
                a_ion.append(a_spec[k])
                ni.append(self.cell2node_bdry(self['ns'][k + 1, :]))
                k_ion.append(k + 1)
            elif spec_list[k] == 'He4':
                print('He4 found')
                nhe4 = self.cell2node_bdry(self['ns'][k + 1, :])
                k_he4 = k + 1
            else:
                print('imp', spec_list[k])
                z_imp.append(z_spec[k])
                a_imp.append(a_spec[k])
                nz.append(self.cell2node_bdry(self['ns'][k + 1, :]))
                k_imp.append(k + 1)

        ni = np.array(ni)
        nz = np.array(nz)
        z_ion = np.array(z_ion)
        a_ion = np.array(a_ion)
        z_imp = np.array(z_imp)
        a_imp = np.array(a_imp)

        n_ion = len(z_ion)
        n_imp = len(z_imp)

        amain = np.sum(np.array([a_ion[k] * ni[k] for k in range(n_ion)]), axis=0)
        zmain = np.sum(np.array([z_ion[k] * ni[k] for k in range(n_ion)]), axis=0)
        nitot = np.sum(np.array([ni[k] for k in range(n_ion)]), axis=0)
        amain /= nitot
        zmain /= nitot

        ne = self.cell2node_bdry(self['ns'][0, :])
        f_imp = np.array([nz[k] / ne for k in range(n_imp)])
        f_ion = np.array([ni[k] / nitot for k in range(n_ion)])

        return {
            'n_ion': n_ion,
            'z_ion': z_ion,
            'a_ion': a_ion,
            'n_imp': n_imp,
            'z_imp': z_imp,
            'a_imp': a_imp,
            'np': ni,
            'nz': nz,
            'nhe4': nhe4,
            'amain': amain,
            'zmain': zmain,
            'f_imp': f_imp,
            'f_ion': f_ion,
            'k_ion': k_ion,
            'k_imp': k_imp,
            'k_he4': k_he4,
            'spec_list': spec_list,
            'z_spec': z_spec,
            'a_spec': a_spec
        }

    def update_particle_balance(self):
        ps_xe = 1.6022e-19
        ps_mp = 1.6726e-27

        z_spec = [round(x) for x in self['qatom_S'][1:] / ps_xe]
        a_spec = [round(x) for x in self['m_S'][1:] / ps_mp]
        n_spec = len(z_spec)

        n_imp = len(self['m_SIMPI'])
        n_ion = n_spec - n_imp

        z_ion = z_spec[0:n_ion]
        a_ion = a_spec[0:n_ion]
        z_imp = z_spec[n_ion:n_imp + n_ion]
        a_imp = a_spec[n_ion:n_imp + n_ion]

        ns_ion = self['ns'][1:n_ion + 1]
        ns_imp = self['ns'][n_ion + 1:n_imp + n_ion + 1]

        nrho = len(self['rho'])
        rho = self['rho'][:]
        ne = self['ns'][0, :]
        density_beam = self['nbeami'][0, :]
        zeff = self['Zeff'][:]

        # ------------------------------------------------------------------
        # density
        density_ion = {}
        for k in range(n_ion):
            density_ion[k] = np.zeros(nrho - 1)

        density_imp = {}
        for k in range(n_imp):
            density_imp[k] = np.zeros(nrho - 1)

        density_th = np.zeros(nrho)

        f_ion = ns_ion / np.sum(ns_ion, axis=0)
        f_imp = ns_imp / np.sum(ns_imp, axis=0)

        for i in range(nrho - 1):
            a = 0
            b = 0
            c = 0
            d = 0
            for k in range(n_imp):
                b = b + f_imp[k, i] * z_imp[k]
                d = d + f_imp[k, i] * z_imp[k] * z_imp[k]
            for k in range(n_ion):
                a = a + f_ion[k, i] * z_ion[k]
                c = c + f_ion[k, i] * z_ion[k] * z_ion[k]

            zne_adj = ne[i]
            zzne_adj = ne[i] * zeff[i]

            # depletion due to beam ions
            zne_adj = zne_adj - 1. * density_beam[i]
            zzne_adj = zzne_adj - 1.**2 * density_beam[i]

            # effective main ion and impurity densities
            nion = (zne_adj * d - zzne_adj * b) / (a * d - b * c)
            nimp = (zzne_adj * a - zne_adj * c) / (a * d - b * c)

            for k in range(n_ion):
                density_ion[k][i] = f_ion[k, i] * nion
            for k in range(n_imp):
                density_imp[k][i] = f_imp[k, i] * nimp

            for k in range(n_ion):
                density_th[i] = density_th[i] + density_ion[k][i]
            for k in range(n_imp):
                density_th[i] = density_th[i] + density_imp[k][i]

        for k in range(n_ion):
            self['ns'][k + 1, :] = density_ion[k]

        for k in range(n_imp):
            self['ns'][k + n_ion + 1, :] = density_imp[k]
