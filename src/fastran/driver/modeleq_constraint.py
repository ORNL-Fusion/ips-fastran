"""
model equilibrium, profile adjust
"""
from fastran.driver.modeleq_constraint_io import update_state, constraint_pedestal_width
from fastran.plasmastate.plasmastate import plasmastate
from fastran.instate import instate_io
from ipsframework import Component

class modeleq_constraint(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print ('Created %s' % (self.__class__))

    def init(self, timestamp=0):
        print ('modeleq_constraint.init() called')

    def step(self, timestamp=0):
        print ('modeleq_constraint.step() called')

        services = self.services
        services.stage_state()

        ps_update = getattr(self, 'PS_UPDATE', 'disabled')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        if ps_update == 'enabled':
            cur_state_file = services.get_config_param('CURRENT_STATE')

        update_state(timestamp,cur_instate_file,nmax_iter=100,const=None)

        k_pedestal_constraint =  int(getattr(self,"PEDESTAL","-1"))

        if k_pedestal_constraint >= 0 and timestamp >= k_pedestal_constraint:
            constraint_pedestal_width(cur_instate_file)

        if ps_update == 'enabled':
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            instate_io.instate_to_ps(cur_instate_file, ps)
            ps.store(cur_state_file)

        services.update_state()
        services.stage_output_files(timestamp, self.OUTPUT_FILES)
