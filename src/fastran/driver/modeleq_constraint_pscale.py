"""
model equilibrium, profile adjust
"""
import numpy as np
from Namelist import Namelist
from ipsframework import Component


class modeleq_constraint_pscale(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('modeleq_constraint_pscale.init() called')
        return

    def step(self, timeid=0):
        print('modeleq_constraint_pscale.step() started')

        self.services.stage_state()
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        pscale = float(self.PSCALE)

        scale_pressure(cur_instate_file, pscale)

        self.services.update_state()
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)


def scale_pressure(f_instate, pscale=1.):
    instate = Namelist(f_instate)
    rho = instate["inmetric"]["rho"]
    pmhd = np.array(instate["instate"]["pmhd"])

    instate["instate"]["pmhd"] = pscale*pmhd

    instate.write(f_instate)
