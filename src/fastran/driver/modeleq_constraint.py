"""
model equilibrium, profile adjust
"""

from ipsframework import Component
from fastran.driver.modeleq_constraint_io import update_state, constraint_pedestal_width
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate


class modeleq_constraint(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('modeleq_constraint.init() called')

    def step(self, timeid=0):
        print('modeleq_constraint.step() started')

        self.services.stage_state()

        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        update_state(cur_instate_file, nmax_iter=100)

        k_pedestal_constraint = int(getattr(self, "PEDESTAL", "-1"))
        if k_pedestal_constraint >= 0 and timeid >= k_pedestal_constraint:
            constraint_pedestal_width(cur_instate_file)

        ps_update = getattr(self, 'PS_UPDATE', 'disabled')
        if ps_update == 'enabled':
            cur_state_file = self.services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
            instate = Instate(cur_instate_file)

            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            instate.to_ps(ps)
            ps.store(cur_state_file)

        self.services.update_state()
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)
