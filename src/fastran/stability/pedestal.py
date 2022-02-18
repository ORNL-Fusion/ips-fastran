"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""

from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate

import numpy as np
from Namelist import Namelist
from fastran.util.modelprofile import profile_pedestal
from fastran.util.zinterp import zinterp
from fastran.util.formula import get_ni
from fastran.util.loglinear import loglinear
from fastran.stability import pedestal_io

class pedestal(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('pedestal.init() called')

    def step(self, timeid=0):
        print('pedestal.step() started')

        #--- stage plasma state files
        self.services.stage_state()

        #--- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #--- eped parameters 
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        fn_inped = getattr(self, 'INPED', 'inped')

        pedestal_io.write_input(cur_instate_file, cur_state_file, ps_backend) 
        pedestal_io.update_state(fn_inped, cur_instate_file, cur_state_file, ps_backend)

        #--- update plasma state files
        self.services.update_state()

        #--- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('pedestal.finalize() called')

