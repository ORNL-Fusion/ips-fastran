"""
 fastran model equilibrium init component
"""

from component import Component
from Namelist import Namelist
from instate_model import instate_model

class modeleq_init(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print(('Created %s' % (self.__class__)))

    def init(self, timestamp=0.0):
        print ('modeleq_init.init() called')

    def step(self, timeStamp):
        print ('modeleq_init.step() called')
        services = self.services

        #-- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #-- plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        #-- update instate from dakota
        fn_instate  = getattr(self, "INSTATE", "instate")
        instate = Namelist(fn_instate)

        var_list = [var.upper() for var in instate["instate"].keys()]
        for key in var_list:
            try:
                instate["instate"][key][0] = float(getattr(self, key))
                print((key, 'updated'))
            except AttributeError:
                pass

        instate.write(cur_instate_file)

        #-- make instate profiles
        instate_model(cur_instate_file)

        #-- touch dummy files
        open(cur_state_file,"w")
        open(cur_eqdsk_file,"w")

        #-- update plasma state
        services.update_plasma_state()

        #-- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0):
        print ('modeleq_init.finalize() called')
