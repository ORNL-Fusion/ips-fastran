"""
 fastran model equilibrium init component
"""

from Namelist import Namelist
from fastran.instate.instate_model import instate_model
from ipsframework import Component


class modeleq_init(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print(('Created %s' % (self.__class__)))

    def init(self, timeid=0):
        print ('modeleq_init.init() called')

    def step(self, timeid=0):
        print ('modeleq_init.step() started')
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
        services.update_state()

        #-- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print ('modeleq_init.finalize() called')
