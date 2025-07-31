"""
 -----------------------------------------------------------------------
 tglfep component
 -----------------------------------------------------------------------
"""

import os
import shutil
from ipsframework import Component
from fastran.stability import tglfep_io


class tglfep(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('tglfep.init() called')

    def step(self, timeid=0):
        print('tglfep.step() started')

        # -- excutable
        tglfep_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(tglfep_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        tglfep_io.write_inputfiles(cur_state_file, cur_eqdsk_file)

        # -- run tglfep
        print('run tglfep')

        cwd = self.services.get_working_dir()
        task_id = self.services.launch_task(self.NPROC, cwd, tglfep_bin, logfile='tglfep.log')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing: tglfep')

        # -- get tglfep output
        tglfep_io.update_state()

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('tglfep.finalize() called')
