#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 genray component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from component import Component
import traceback
#-------------------
#--- zcode libraries
import zgenray

class genray(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        services = self.services
        services.stage_plasma_state()

        pstool_path = self.PSTOOL_PATH

        pstool_bin = pstool_path

        print 'pstool_bin path:'
        print pstool_bin

        print 'pstool genray init'

        logfile=open("pstool_init.log","w")
        retcode = subprocess.call([pstool_bin, "init", "genray"],
                      stdout=logfile,stderr=logfile)
        logfile.close()
        if (retcode != 0):
            logMsg = 'Error executing ' + pstool_bin
            services.exception(logMsg)
            raise Exception(logMsg)

        cur_state_file = services.get_config_param('CURRENT_STATE')
        try:
            shutil.copyfile("ips-state-genray.nc",cur_state_file)
        except Exception:
            logMsg = 'ERROR : File copy failed'
            services.exception(logMsg)
            raise

        try:
            services.update_plasma_state()
        except Exception:
            logMsg = 'Error in call to update_plasma_state()'
            services.exception(logMsg)
            raise

        return

    def step(self, timeStamp=0.0):

        #--- entry

        if (self.services == None) :
            logMsg = 'Error in genray.step() : No services'
            raise Exception(logMsg)
        services = self.services

        #--- excutable

        try:
            genray_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            genray_bin = os.path.join(self.BIN_PATH, 'xgenray')
            logMsg = 'Using default binary path for GENRAY'
            services.info(logMsg)
        print genray_bin

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
        except Exception:
            logMsg = 'Error in call to stage_plasma_state()'
            services.exception(logMsg)
            raise

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- generate genray input

        f_ingenray = "ingenray"
        zgenray.io_write_inputfiles(cur_state_file,cur_eqdsk_file,f_ingenray)

        #--- run genray

        print 'run genray'

        #logfile=open("genray.log","w")
        #retcode = subprocess.call([genray_bin],
        #              stdout=logfile,stderr=logfile)
        #logfile.close()

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, genray_bin, logfile = 'xgenray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           logMsg = 'Error executing ', 'xgenray'
           services.exception(logMsg)
           raise Exception(logMsg)

        #--- get genray output

        zgenray.io_update_state(cur_state_file,cur_eqdsk_file)

        #--- update plasma state files

        try:
            services.update_plasma_state()
        except Exception, e:
            logMsg =  'Error in call to update_plasma_state()'
            services.exception(logMsg)
            raise

        #--- archive output files

        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception: 
            logMsg = 'Error in call to stage_output_files()'
            services.exception(logMsg)
            raise Exception

        return


    def restart(self, timeStamp=0.0):

        self.services.info('restart() called for GENRAY')
        restart_root = self.services.get_config_param('RESTART_ROOT')
        self.services.get_restart_files(restart_root, timeStamp, self.RESTART_FILES)

        return
    

    def checkpoint(self, timestamp=0.0):

        services = self.services

        services.info('checkpoint() called for genray')

        services.save_restart_files(timestamp, self.RESTART_FILES)

        return

    def finalize(self, timeStamp=0.0):
        return
    
