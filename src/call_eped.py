# nested workflow wrapper for EPED

from  component import Component

#-------------------
#--- zcode libraries
from Namelist import Namelist

class call_eped(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        print 'entry: call_eped init'

        (sim_name, init, driver) = \
            self.services.create_sub_workflow('test_flow', self.SUB_WORKFLOW, {})
        self.init = init
        self.driver = driver
        print sim_name, init, driver

        self.services.stage_input_files(self.INPUT_FILES)

        self.services.call(self.init, 'step', timeStamp)
        self.services.call(self.driver,'init', timeStamp)

        print "success***"

        return

    def validate(self, timeStamp=0):

        return

    def step(self, timeStamp=0):

        print 'entry: call_eped'

        #--- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')

        #--- stage plasma state files
        self.services.stage_plasma_state()

        #--- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #--- componse eped input file form the plasma state
        io_input_eped(cur_state_file)

        #--- call nested eped workflow, driver step
        self.services.call(self.driver, 'step', '0.0')

        #--- stage output files (output files from sub workflow driver)
        self.services.stage_subflow_output_files()

        #--- update plasma state from eped output
        f_eped_output = 'eped_state.nc' # need to access subwork flow configuration variable...
        io_update_state(cur_state_file,f_eped_output)

        return

    def finalize(self, timeStamp=0.0):

        self.services.call(self.driver, 'finalize', '0.0')

        return

def io_input_eped(f_fstate):

    print 'entry: i0_input_eped'

    # generate eped.in from plasma state_file
    # more to come...

def io_update_state(cur_state_file,f_eped_output):

    print 'entry: i0_update_state()'

    # update plasma state from eped output
    # more to come...
