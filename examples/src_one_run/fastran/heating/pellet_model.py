"""
 -----------------------------------------------------------------------
 model pellet component
 -----------------------------------------------------------------------
"""

import numpy as np
from scipy import optimize
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.util import dakota_io
from fastran.util.fastranutil import namelist_default


class pellet_model(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('pellet_model.init() called')

    def step(self, timeid=0):
        print('pellet_model.init() started')

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- idakota binding
        inpellet = Namelist("inpellet", case="upper")

        print("start dakota update")
        dakota_io.update_namelist(self, inpellet)
        inpellet.write("inpellet_updated")

        # -- pellet model profile
        self.update_state(cur_state_file, inpellet)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('pellet_model.finalize() called')

    def gauss(self, p, x):
        return p[0]*np.exp(-(x - p[1])**2/(2*p[2]**2))

    def gauss_profile(self, rho, vol, xmid, xwid):
        nrho = len(rho)
        def func(y0):
            y = self.gauss([y0, xmid, xwid], rho)
            integral = 0.
            for i in range(nrho - 1):
                integral += 0.5*(y[i+1] + y[i])*(vol[i + 1] - vol[i])
            return integral - 1.
        f0 = optimize.bisect(func, 0., 10., xtol=1.e-6)
        profile = self.gauss([f0, xmid, xwid], rho)
        return profile

    def update_state(self, cur_state_file, inpellet):
        npellet = namelist_default(inpellet, "inpellet", "npellet", [0])[0]
        sn_peak = namelist_default(inpellet, "inpellet", "sn_peak", [0.])
        sn_total = namelist_default(inpellet, "inpellet", "sn_total", [0.])
        xmid = namelist_default(inpellet, "inpellet", "xmid", [0.])
        xwid = namelist_default(inpellet, "inpellet", "xwid", [0.1])

        print("npellet =", npellet)
        print("sn_peak =", sn_peak)
        print("xmid =", xmid)
        print("xwid =", xwid)

        # -- read ps
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        rho = ps["rho"][:]
        vol = ps["vol"][:]
        nrho = len(rho)

        # -- profile
        sn_sum = np.zeros(nrho)
        for k in range(npellet):
            profile = self.gauss_profile(rho, vol, xmid[k], xwid[k])
            print(profile)
            if sn_total[k] > 0.:
               sn_sum = sn_sum + profile*sn_total[k]*1.e19
            else:
               sn_sum = sn_sum + profile*sn_peak[k]*1.e19
            print(sn_sum)
        # ps.load_vol_profile(rho, sn_sum, "rho_nbi", "sbsce", k=0)
        ps.load_vol_profile(rho, sn_sum, "rho", "sn_trans", k=0) # note that sn_trans is not consistent with the definition of plasma state

        # -- write to ps
        ps.store(cur_state_file)
