#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 curray component 
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from  component import Component

#--- zcode libraries
import zcurray

class curray(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        services = self.services

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- excutable

        try:
            curray_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            curray_bin = os.path.join(self.BIN_PATH, 'xcurray')

        #--- generate curray input

        zcurray.wrt_curray_input(cur_state_file,"incurray",cur_eqdsk_file)

        #--- run curray

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, curray_bin, "curray_in", logfile='xcurray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'curray'
           raise

        #--- get curray output

        zcurray.read_curray_output()

        zcurray.io_update_state(cur_state_file,cur_eqdsk_file,"outcurray","incurray")

        #--- update plasma state files

        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            raise

        #--- archive output files

        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            raise Exception, e

        return
    
    def finalize(self, timeStamp=0.0):
        return
    
