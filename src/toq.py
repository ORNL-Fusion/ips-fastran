"""
 -----------------------------------------------------------------------
 toq component
 -----------------------------------------------------------------------
"""

import os
import shutil
from component import Component
from Namelist import Namelist
from plasmastate import plasmastate

class toq(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0.0):
        print('toq.init() called')

    def step(self, timeStamp=0.0):
        print('enter toq.step()')

        #--- code entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #-- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- write input
        f_intoq = getattr(self,'INTOQ','intoq')

        try:
            shutil.copy(f_intoq, 'intoq')
        except shutil.SameFileError:
            pass

        intoq = Namelist(f_intoq)
        intoq["input"]["fneqdsk"] = [cur_eqdsk_file]
        intoq.write(f_intoq)

        #--- excutables
        toq_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(toq_bin)

        #--- run toq
        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, toq_bin, logfile = 'xtoq.log')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing toq')

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0.0):
        print('enter toq.step()')
