#!/usr/bin/env python

"""
model equilibrium, profile adjust
"""

from component import Component

from modeleq_constraint_io import init_plasmastate
from modeleq_constraint_io import update_state, constraint_pedestal_width

class modeleq_constraint(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):

        print 'modeleq_constraint.init() called'

        services = self.services
        services.stage_plasma_state()
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        init_plasmastate(cur_instate_file)

        services.update_plasma_state()

        return

    def step(self, timestamp=0):

        #--- entry

        services = self.services
        services.stage_plasma_state()
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        update_state(timestamp,cur_instate_file,nmax_iter=100,const=None)

        k_pedestal_constraint =  int(getattr(self,"PEDESTAL","-1"))

        if k_pedestal_constraint >= 0 and timestamp >= k_pedestal_constraint:
            constraint_pedestal_width(cur_instate_file)

        services.update_plasma_state()
        services.stage_output_files(timestamp, self.OUTPUT_FILES)
