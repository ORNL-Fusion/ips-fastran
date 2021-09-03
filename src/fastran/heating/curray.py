"""
 -----------------------------------------------------------------------
 curray component
 -----------------------------------------------------------------------
"""

import os
from fastran.heating import curray_io
from ipsframework import Component

class curray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('curray.init() called')

    def step(self, timeid=0):
        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- excutable
        curray_bin = os.path.join(self.BIN_PATH, self.BIN)

        #--- generate curray input
        curray_io.wrt_curray_input(cur_state_file, "incurray", cur_eqdsk_file)

        #--- run curray
        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, curray_bin, "curray_in", logfile='xcurray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           raise Exception('Error executing curray')

        #--- get curray output
        curray_io.read_curray_output()
        curray_io.update_state(cur_state_file, cur_eqdsk_file, "outcurray", "incurray")

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0.0):
        print('curray.finalize() called')
