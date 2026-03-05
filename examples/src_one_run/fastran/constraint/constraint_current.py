"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""
from fastran.plasmastate.plasmastate import plasmastate
from fastran.constraint import constraint_current_io
from ipsframework import Component
from fastran.util.fastranutil import freeze

class constraint_current(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('constraint_current.init() done')

    def step(self, timeid=0):
        print('constraint_current.step() started')

        # -- freeze/resume
        if freeze(self, timeid, 'constraint_current'): return None

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- apply constraint to local plasma state
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
            jboot = float(getattr(self, 'JBOOT', '1.0'))

            print('q0 = ', q0)
            print('alpha = ', alpha)

            constraint_current_io.parabolicJ(cur_state_file, cur_eqdsk_file, rho_jbdry=rho_jbdry, alpha=alpha, q0=q0, jboot=jboot)

        elif self.METHOD == 'hat':
            q0 = float(self.Q0)
            jboot = float(getattr(self, 'JBOOT', '1.0'))

            wid_jhat = float(getattr(self, 'WID_JHAT', '0.2'))
            rho_jhat_min = float(getattr(self, 'RHO_JHAT_MIN', '0.1'))
            rho_jhat_max = float(getattr(self, 'RHO_JHAT_MAX', '0.6'))

            print('****', wid_jhat)
            print('****', rho_jhat_min)
            print('****', rho_jhat_max)

            constraint_current_io.hatJ(
                cur_state_file,
                cur_eqdsk_file,
                rho_jbdry=rho_jbdry,
                q0=q0,
                wid_jhat=wid_jhat,
                rho_jhat_min=rho_jhat_min,
                rho_jhat_max=rho_jhat_max,
                niter_max=100,
                jboot=jboot
            )

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('constraint_current.finalize() called')
