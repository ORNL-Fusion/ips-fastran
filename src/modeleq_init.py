#! /usr/bin/env python

"""
 fastran model equilibrium init component
"""

from component import Component
from Namelist import Namelist

class modeleq_init(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('modeleq_init.init() called')

    def step(self, timeStamp):

        #-- entry

        print ('modeleq_init.step() called')
        services = self.services

        #-- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #-- plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #-- update instate from dakota

        fn_instate  = getattr(self,"INSTATE","instate")

        instate = Namelist(fn_instate)

        var_list = [ var.upper() for var in instate["instate"].keys()]
        for key in var_list:
            try:
                instate["instate"][key][0] = float(getattr(self, key))
                print key, 'updated'
            except AttributeError:
                pass

        instate.write(cur_instate_file)

        inbc = Namelist()
        inbc["inbc"]["r0"] = instate["instate"]["r0"]
        inbc["inbc"]["b0"] = instate["instate"]["b0"]
        inbc["inbc"]["ip"] = instate["instate"]["ip"]
        inbc["inbc"]["nbdry"] = instate["instate"]["nbdry"]
        inbc["inbc"]["rbdry"] = instate["instate"]["rbdry"]
        inbc["inbc"]["zbdry"] = instate["instate"]["zbdry"]
        inbc.write(cur_bc_file)

        #-- touch dummy files

        open(cur_state_file,"w")
        open(cur_eqdsk_file,"w")
        open(cur_fastran_file,"w")

        #-- update plasma state

        services.update_plasma_state()

        #-- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0):

        print 'modeleq_init.finalize() called'
