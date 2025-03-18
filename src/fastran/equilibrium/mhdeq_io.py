'''
 -----------------------------------------------------------------------
  mhdeq_io
 -----------------------------------------------------------------------
'''

import shutil
import numpy as np
from scipy.interpolate import interp1d
from Namelist import Namelist
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver.inmetric_io import ps_to_inmetric
from fastran.state.instate import Instate
from fastran.util.zinterp import zinterp


class mhdeq_io():
    def __init__(self, model_pressure):
        self.model_pressure = model_pressure
        pass

    def from_instate(self, instate=None, f_instate=''):
        if instate:
            self.instate = instate
        elif f_instate:
            self.instate = Instate(f_instate)
        else:
            self.instate = None

    def from_ps(self, f_ps):
        self.instate.from_ps(f_ps)

    def get_p(self):
        if self.model_pressure == 'kinetic':
            print('model_pressure = kinetic')
            ne = self.instate['ne']
            te = self.instate['te']
            ti = self.instate['ti']
            wbeam = self.instate['wbeam']
            walp = self.instate['walpha']
            zeff = self.instate['zeff']
            nth = self.instate['density_th']
            pmhd = 1.602e3 * (ne * te + nth * ti) + 2. / 3. * 1.e6 * (wbeam + walp)
        elif self.model_pressure == 'totol':
            print('model_pressure = total')
            pmhd = instate['p_eq']
        else:
            raise Exception(f'input error: model_pressure = {model_pressure} not supported')
        return pmhd

    def get_jpar(self):
        return self.instate['j_tot']

    def load_geqdsk_relax(self, f_geqdsk, relax=0):
        if relax > 0:
            inmetric0 = self.inmetric
        self.load_geqdsk(f_geqdsk)
        if relax > 0:
            print('relax metric', relax)
            for key in self.inmetric:
                self.inmetric[key] = (1. - relax) * inmetric0[key] + relax * self.inmetric[key]

    def load_geqdsk(self, f_geqdsk):
        # nrho for ps grid of geqdsk, i.e. corresponding to number of contours
        # used inmetric for jpar => ffprim: psi, rho, rhob, ipol, vol, gr2
        nrho = self.instate['nrho'][0]
        self.instate.inmetric(f_geqdsk, nrho=nrho, nth=101)
        self.inmetric = {key: self.instate.data['inmetric'][key]
                         for key in ['nrho', 'psi', 'rho', 'rhob', 'ipol', 'vol', 'gr2i']}

    def get_ffprime_pprime_relax(self, relax=0, max_error=0.01):
        if relax > 0:
            pprim0  = self.pprim
            ffprim0 = self.ffprim
            f0      = self.f
            p0      = self.p
        self.get_ffprime_pprime()
        iconv = 0
        if relax > 0:
            self.pprim  = (1. - relax) * pprim0  + relax * self.pprim
            self.ffprim = (1. - relax) * ffprim0 + relax * self.ffprim
            self.f      = (1. - relax) * f0      + relax * self.f
            self.p      = (1. - relax) * p0      + relax * self.p

            diff_p = np.sum(self.pprim - pprim0) / np.sum(self.pprim)
            diff_f = np.sum(self.ffprim - ffprim0) / np.sum(self.ffprim) 
            print (f'diff_p, diff_f = {diff_p}, {diff_f}')
            if diff_p < max_error and diff_f < max_error:
                iconv = 1
        return iconv
           
    def get_ffprime_pprime(self):

        mu0 = 4.e-7 * np.pi
        r0 = self.instate['r0'][0]  # [m]
        b0 = self.instate['b0'][0]  # [T]
        ip = self.instate['ip'][0] * 1.e6  # [A]
        b0 = abs(b0)

        rho_in = self.instate['rho']  # []
        press_in = self.get_p()  # [Pa]
        jpar_in = self.get_jpar() * 1.e6  # [A/m^2]

        psirho = self.inmetric['psi'] / self.inmetric['psi'][-1]
        dpsi = abs(self.inmetric['psi'][0] - self.inmetric['psi'][-1])

        rhob = self.inmetric['rhob'][0]
        nrho = self.inmetric['nrho'][0]
        rho = self.inmetric['rho']
        drho = 1. / (nrho - 1.)

        npsi = nrho
        psi = np.arange(npsi) / (npsi - 1.)

        rhopsi_spl = interp1d(psirho, rho, kind='cubic')
        rho_eval = rhopsi_spl(psi)
        rho_eval[0] = 0.
        rho_eval[-1] = 1.

        ipol = self.inmetric['ipol']
        vol = self.inmetric['vol']
        rm2 = self.inmetric['gr2i']

        rm2_spl = zinterp(rho, rm2)
        rm2_psi = rm2_spl(rho_eval)

        # -- jpar to jtor
        jpar = zinterp(rho_in, jpar_in)(rho)

        curt = np.zeros(nrho)
        curt[0] = 0.
        for i in range(1, nrho):
            dV = vol[i] - vol[i - 1]
            ipol_m = 0.5 * (ipol[i] + ipol[i - 1])
            jpar_m = 0.5 * (jpar[i] + jpar[i - 1])
            curt[i] = (curt[i - 1] / ipol[i - 1] + jpar_m / ipol_m**2 / (2. * np.pi * r0) * dV) * ipol[i]

        jtor = np.zeros(nrho + 1)
        rho_cell = np.zeros(nrho + 1)
        for i in range(1, nrho):
            dV = vol[i] - vol[i - 1]
            jtor[i] = (curt[i] - curt[i - 1]) / dV * 2. * np.pi * r0
            rho_cell[i] = 0.5 * (rho[i] + rho[i - 1])
        jtor[0] = jtor[1]
        jtor[-1] = 2. * jtor[-2] - jtor[-3]
        rho_cell[0] = 0.
        rho_cell[-1] = 1.

        jtor_spl = zinterp(rho_cell, jtor)
        jtor_psi = jtor_spl(rho_eval)

        press_spl = zinterp(rho_in, press_in, s=0)
        press_psi = press_spl(rho_eval, der=0)

        pprim_psi = zinterp(psi, press_psi, s=0)(psi, der=1)
        pprim_psi = pprim_psi / dpsi

        ipol_spl = zinterp(rho, ipol)
        ipol_psi = ipol_spl(rho_eval, der=0)  # ------

        # -- Ip scale
        ip_temp = 0.0
        for i in range(nrho - 1):
            dV = vol[i + 1] - vol[i]
            rho_m = 0.5 * (rho[i + 1] + rho[i])
            ip_temp += dV * jtor_spl(rho_m)
        ip_temp /= 2. * np.pi * r0
        print('ip, ip cal, ip/(ip cal) = %5.3f %5.3f %5.3f' % (ip * 1.e-6, ip_temp * 1.e-6, ip / ip_temp))

        for i in range(len(jtor_psi)):
            jtor_psi[i] = ip / ip_temp * jtor_psi[i]

        # -- calculate p' and ff'
        x0 = [0.1, 0.15, 0.2]
        extrapolation(psi, pprim_psi, x0, y0=None)

        ffprim_psi = np.zeros(npsi)
        for i in range(npsi):
            ffprim_psi[i] = -mu0 * r0 / (r0**2 * rm2_psi[i]) * (jtor_psi[i] + r0 * pprim_psi[i])

        pprim_psi[-1] = extrapolation_sep(
            [psi[-4], psi[-3], psi[-2]],
            [pprim_psi[-4], pprim_psi[-3], pprim_psi[-2]])

        ffprim_psi[-1] = extrapolation_sep(
            [psi[-4], psi[-3], psi[-2]],
            [ffprim_psi[-4], ffprim_psi[-3], ffprim_psi[-2]])

        # -- output
        self.npsi = [npsi]
        self.psin = psi
        self.pprim = -pprim_psi
        self.ffprim = -ffprim_psi
        self.f = ipol_psi * r0 * b0
        self.p = press_psi


# -- extrapolation routines

def inverse3x3(x, y):
    a = np.array(
        [[x[0]**2, x[0], 1.],
         [x[1]**2, x[1], 1.],
         [x[2]**2, x[2], 1.]])
    a_inv = np.linalg.inv(a)
    return np.dot(a_inv, y)


def cubic(p, x0):
    return p[0] * x0**2 + p[1] * x0 + p[2]


def extrapolation(rho, f, x0, y0):
    spl = zinterp(rho, f, s=0)

    if y0 is None:
        # p = inverse3x3(x0, [spl(x) for x in x0])
        tmp = spl(x0, der=0)
        p = inverse3x3(x0, tmp)
        for k in range(len(rho)):
            if rho[k] < x0[0]:
                f[k] = cubic(p, rho[k])
    else:
        # p = inverse3x3(x0, [y0, spl(x0[1]), spl(x0[2])])
        tmp = spl([x0[1], x0[2]], der=0)
        p = inverse3x3(x0, [y0, tmp[0], tmp[1]])
        for k in range(len(rho)):
            if rho[k] < x0[1]:
                f[k] = cubic(p, rho[k])


def extrapolation_sep(x0, y0):
    p = inverse3x3(x0, y0)
    return cubic(p, 1.)
