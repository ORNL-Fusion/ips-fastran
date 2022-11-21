"""
 -----------------------------------------------------------------------
 model h/cd component
 -----------------------------------------------------------------------
"""

import numpy as np
from scipy import optimize
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.equilibrium.efit_eqdsk import readg
from fastran.util import dakota_io


class hcd_model(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('hcd_model.init() called')

    def step(self, timeid=0):
        print('hcd_model.init() started')

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        if ps_backend == 'INSTATE':
            cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- idakota binding
        inhcd = Namelist("inhcd", case="upper")

        print("start dakota update")
        dakota_io.update_namelist(self, inhcd)
        inhcd.write("inhcd_updated")

        add = int(getattr(self, "ADD", "0"))
        print('add = ', add)
        print('ps_backend = ', ps_backend)

        if ps_backend == 'PS':
            self.update_state(cur_state_file, cur_eqdsk_file, inhcd, add)
        elif ps_backend == 'INSTATE':
            self.update_instate(cur_instate_file, cur_eqdsk_file, inhcd, add)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('hcd_model.finalize() called')

    def gauss_profile(self, rho, vol, xmid, xwid):
        nrho = len(rho)
        def gauss(p, x): return p[0]*np.exp(-(x - p[1])**2/(2*p[2]**2))

        def func(y0):
            y = gauss([y0, xmid, xwid], rho)
            P = 0.
            for i in range(nrho - 1):
                P += 0.5*(y[i+1] + y[i])*(vol[i + 1] - vol[i])
            return P - 1.0
        f0 = optimize.bisect(func, 0., 10., xtol=1.e-6)
        profile = gauss([f0, xmid, xwid], rho)
        return profile

    def update_state(self, cur_state_file, cur_eqdsk_file, inhcd, add):
        nsrc = inhcd["inhcd"]["nsrc"][0]
        Pe = inhcd["inhcd"]["Pe"]
        Pi = inhcd["inhcd"]["Pi"]
        xmid = inhcd["inhcd"]["xmid"]
        xwid = inhcd["inhcd"]["xwid"]

        j0_seed = inhcd["inhcd"]["j0_seed"]
        x0_seed = inhcd["inhcd"]["x0_seed"]
        drho_seed = inhcd["inhcd"]["drho_seed"]

        print("Pe =", Pe)
        print("Pi =", Pi)
        print("xmid =", xmid)
        print("xwid =", xwid)

        # -- read ps
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        geq = readg(cur_eqdsk_file)
        r0 = geq["rzero"]
        b0 = abs(geq["bcentr"])
        ip = geq['cpasma']

        rho = ps["rho"][:]
        vol = ps["vol"][:]
        nrho = len(rho)

        #-- heating
        def vol_integral(f):
            val = 0.
            for i in range(nrho - 1):
                f_m = 0.5*(f[i + 1] + f[i])
                dvol = vol[i + 1] - vol[i]
                val += f_m*dvol
            return val

        pe_sum = np.zeros(nrho)
        pi_sum = np.zeros(nrho)
        for k in range(nsrc):
            profile = self.gauss_profile(rho, vol, xmid[k], xwid[k])
            pe_sum = pe_sum + profile*Pe[k]*1.e6
            pi_sum = pi_sum + profile*Pi[k]*1.e6

        ps.load_vol_profile(rho, pe_sum, "rho_icrf", "picrf_totals", k=0, add=add)
        ps.load_vol_profile(rho, pi_sum, "rho_icrf", "picrf_totals", k=1, add=add)

        #-- current
        j_sum = np.zeros(nrho)
        for k in range(nsrc):
            if j0_seed[k] > 0.0:
                j_seed = j0_seed[k]*np.exp(-(rho - x0_seed[k])**2/(2*drho_seed[k]**2))
            else:
                j_seed = np.zeros(nrho)
            j_sum += j_seed
        ps.load_j_parallel(rho, j_sum*1.0e6, "rho_icrf", "curich", r0, b0, add=add)

        # -- write to ps
        ps.store(cur_state_file)

    def update_instate(self, cur_instate_file, cur_eqdsk_file, inhcd, add):
        nsrc = inhcd["inhcd"]["nsrc"][0]
        Pe = inhcd["inhcd"]["Pe"]
        Pi = inhcd["inhcd"]["Pi"]
        xmid = inhcd["inhcd"]["xmid"]
        xwid = inhcd["inhcd"]["xwid"]

        j0_seed = inhcd["inhcd"]["j0_seed"]
        x0_seed = inhcd["inhcd"]["x0_seed"]
        drho_seed = inhcd["inhcd"]["drho_seed"]

        print("Pe =", Pe)
        print("Pi =", Pi)
        print("xmid =", xmid)
        print("xwid =", xwid)

        # -- read instate
        instate = Namelist(cur_instate_file)

        rho = np.array(instate["inmetric"]["rho"])
        rhob = instate["inmetric"]["rhob"][0]
        volp = instate["inmetric"]["volp"]
        nrho = len(rho)

        vol = np.zeros(nrho)
        for i in range(nrho - 1):
            vol[i + 1] = vol[i] + 0.5*(volp[i] + volp[i + 1])*(rho[i + 1] - rho[i])*rhob

        def vol_integral(f):
            val = 0.
            for i in range(nrho - 1):
                f_m = 0.5*(f[i + 1] + f[i])
                dvol = vol[i + 1] - vol[i]
                val += f_m*dvol
            return val

        #-- heating
        pe_sum = np.zeros(nrho)
        pi_sum = np.zeros(nrho)
        for k in range(nsrc):
            profile = self.gauss_profile(rho, vol, xmid[k], xwid[k])
            pe_sum = pe_sum + profile*Pe[k]
            pi_sum = pi_sum + profile*Pi[k]

        if add:
            instate["instate"]["pe_ic"] = np.array(instate["instate"]["pe_ic"]) + pe_sum
            instate["instate"]["pi_ic"] = np.array(instate["instate"]["pi_ic"]) + pi_sum
        else:
            instate["instate"]["pe_ic"] = pe_sum
            instate["instate"]["pi_ic"] = pi_sum

        #-- current
        j_sum = np.zeros(nrho)
        for k in range(nsrc):
            if j0_seed[k] > 0.:
                j_seed = j0_seed[k]*np.exp(-(rho - x0_seed[k])**2/(2*drho_seed[k]**2))
            else:
                j_seed = np.zeros(nrho)
            j_sum += j_seed
        instate["instate"]["j_ic"] = j_sum

        # -- write to instate
        instate.write(cur_instate_file)
