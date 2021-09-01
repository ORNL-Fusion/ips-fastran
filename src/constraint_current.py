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
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('constraint_current.init() done')

    def step(self, timeid=0):
        print('constraint_current.step() started')

        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- freeze
        start_time = int(getattr(self, "START_TIME", "-1"))
        end_time = int(getattr(self, "END_TIME", "1000"))
        if timeid < start_time: return
        if timeid > end_time: return

        #--- apply constraint to local plasma state
        rho_jbdry = float(getattr(self, "RHO_JBDRY", "0.85"))

        if self.METHOD == 'broadJ':
            jaxis = float(self.JAXIS)
            rho_j = float(self.RHO_J)

            print('jaxis =', jaxis)
            print('rho_j =', rho_j)

            constraint_current_io.broadJ_spline(cur_state_file, cur_eqdsk_file, rho_j, jaxis, rho_jbdry)

        elif self.METHOD == 'broadJ_Gauss':
            rho_j = float(self.RHO_J)
            j_in = float(self.J_IN)
            j_out = float(self.J_OUT)

            print('rho_j = ', rho_j)
            print('j_in= ', j_in)
            print('j_out= ', j_out)

            constraint_current_io.broadJ(cur_state_file, cur_eqdsk_file, rho_j, j_in, j_out, rho_jbdry)

        elif self.METHOD == 'parabolicJ':
            q0 = float(self.Q0)
            alpha = float(self.ALPHA)

            print('q0 = ', q0)
            print('alpha = ', alpha)

            constraint_current_io.parabolicJ(cur_state_file, cur_eqdsk_file, rho_jbdry=rho_jbdry, alpha=alpha, q0=q0)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('constraint_current.finalize() called')
