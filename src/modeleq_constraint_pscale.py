"""
model equilibrium, profile adjust
"""

from component import Component
from Namelist import Namelist
from numpy import *

class modeleq_constraint_pscale(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('modeleq_constraint_pscale.init() called')
        return

    def step(self, timeid=0):
        print('modeleq_constraint_pscale.step() started')

        #--- entry
        services = self.services
        services.stage_plasma_state()
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        pscale = float(self.PSCALE)

        scale_pressure(cur_instate_file, pscale)

        services.update_plasma_state()
        services.stage_output_files(timeid, self.OUTPUT_FILES)

def scale_pressure(f_instate, pscale=1.0):

    instate = Namelist(f_instate)
    rho = instate["inmetric"]["rho"]
    pmhd = array(instate["instate"]["pmhd"])

    instate["instate"]["pmhd"] = pscale*pmhd

    instate.write(f_instate)
