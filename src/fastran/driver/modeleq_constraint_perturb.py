"""
model equilibrium, profile adjust
"""
import numpy as np
from scipy.optimize import bisect
from Namelist import Namelist
from fastran.util.zinterp import zinterp
from ipsframework import Component


class modeleq_constraint_perturb(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('modeleq_constraint_perturb.init() called')

    def step(self, timeid=0):
        print('modeleq_constraint_perturb.step() started')

        self.services.stage_state()
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        width = float(self.WIDTH)
        jpert = float(self.JPERT)
        ppert = float(self.PPERT)

        perturb(cur_instate_file, width, jpert, ppert)

        self.services.update_state()
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)


def perturb(f_instate, w=0.3, jpert=0., ppert=0.):
    instate = Namelist(f_instate)
    rho = instate["inmetric"]["rho"]
    qmhd = instate["inmetric"]["qmhd"]
    j_tot = np.array(instate["instate"]["j_tot"])
    pmhd = np.array(instate["instate"]["pmhd"])
    nrho = len(rho)

    qmhd_s = zinterp(rho, qmhd)
    def q2(x):
        return qmhd_s(x) - 2.0

    rho_min = rho[np.argmin(qmhd)]
    rho_2 = bisect(q2, rho_min, 0.8)
    print(rho_2, rho_min, qmhd_s(rho_2))

    j_tot_2 = zinterp(rho, j_tot)(rho_2)
    pmhd_2 = zinterp(rho, pmhd)(rho_2)

    jadd = np.zeros(nrho)
    padd = np.zeros(nrho)
    for i in range(nrho):
        if rho[i] < rho_2-w or rho[i] > rho_2+w :
            jadd[i] = 0.0
            padd[i] = 0.0
        else:
            jadd[i] = sin( (rho[i]-rho_2)/w*pi )
            padd[i] = sin( (rho[i]-rho_2)/w*pi )
    j_tot = j_tot + jadd*jpert*j_tot_2
    pmhd = pmhd + padd*ppert*pmhd_2

    instate["instate"]["j_tot"] = j_tot
    instate["instate"]["pmhd"] = pmhd

    instate.write(f_instate)
