#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 curray component 
 -----------------------------------------------------------------------
"""

import os
from  component import Component
import curray_io

class curray(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        return

    def step(self, timeStamp=0):

        #--- entry

        services = self.services

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- excutable

        curray_bin = os.path.join(self.BIN_PATH, self.BIN)

        #--- generate curray input

        curray_io.wrt_curray_input(cur_state_file,"incurray",cur_eqdsk_file)

        #--- run curray

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, curray_bin, "curray_in", logfile='xcurray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           raise Exception('Error executing curray')

        #--- get curray output

        curray_io.read_curray_output()

        curray_io.update_state(cur_state_file,cur_eqdsk_file,"outcurray","incurray")

        #--- update plasma state files

        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return
    
    def finalize(self, timeStamp=0.0):
        return
    
