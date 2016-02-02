#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 genray component 
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from component import Component

#--- zcode libraries
import zgenray

class genray(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        services = self.services

        #--- excutable

        try:
            genray_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            genray_bin = os.path.join(self.BIN_PATH, 'xgenray')
        print genray_bin

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

        #--- generate genray input

        f_ingenray = "ingenray"
        zgenray.io_write_inputfiles(cur_state_file, cur_eqdsk_file, f_ingenray)

        #--- run genray

        print 'run genray'

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, genray_bin, logfile = 'xgenray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'xgenray'
           raise

        #--- get genray output

        zgenray.io_update_state(cur_state_file, cur_eqdsk_file)

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
    
