#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""

from component import Component
from plasmastate import plasmastate

import constraint_current_io

class constraint_current(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        services = self.services

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file name

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- apply constraint to local plasma state

        rho_jbdry = float(getattr(self,"RHO_JBDRY","0.85"))

        if self.METHOD == 'broadJ':

            # rho_j = float(self.RHO_J)
            # j_in = float(self.J_IN)
            # j_out = float(self.J_OUT)
            #
            # print 'rho_j = ', rho_j
            # print 'j_in= ', j_in
            # print 'j_out= ', j_out
            #
            # constraint_current_io.broadJ(cur_state_file, cur_eqdsk_file, rho_j, j_in, j_out,rho_jbdry)

            jaxis = float(self.JAXIS)
            rho_j = float(self.RHO_J)

            constraint_current_io.broadJ_spline(cur_state_file, cur_eqdsk_file, rho_j, jaxis, rho_jbdry)

        elif self.METHOD == 'parabolicJ':

            q0 = float(self.Q0)
            alpha = float(self.ALPHA)

            print 'q0 = ', q0
            print 'alpha = ', alpha

            constraint_current_io.parabolicJ(cur_state_file, cur_eqdsk_file, cur_bc_file, rho_jbdry = rho_jbdry, alpha=alpha, q0=q0)

        #--- update plasma state files

        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def finalize(self, timeStamp=0.0):

        return