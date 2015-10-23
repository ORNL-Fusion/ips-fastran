#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran transport solver component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from  component import Component 

#-------------------
#--- zcode libraries
import zfastran

class fastran_solver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        if (self.services == None) :
            print 'Error in fastran_solver.step() : No services'
            raise Exception('Error in fastran_solver.step(): No services')
        services = self.services

        #--- stage plasma state files 

        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e

        #--- get plasma state file name

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        #cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        try:
            shutil.copyfile(self.INFASTRAN,"infastran")
        except:
            pass

        #--- generate fastran input

        zfastran.io_write_input(cur_state_file,cur_eqdsk_file)

        #--- run fastran

        try:
            fastran_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            fastran_bin = os.path.join(self.BIN_PATH, 'xfastran')
        print fastran_bin

        ncpu =  int(self.NPROC)
        nky  =  int(self.NPROC_KY)
        n1d  =  ncpu/nky

        print "ncpu = ",ncpu
        print "n1d  = ",n1d 
        print "nky  = ",nky

        cwd = services.get_working_dir()
        task_id = services.launch_task(ncpu, cwd, fastran_bin, "%d"%n1d, "%d"%nky, logfile = 'xfastran.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'fastran'
           raise

        #--- update local plasma state
        try:
           relax = float(self.RELAX)
        except:
           relax = 0.5
        print 'relax=',relax

        zfastran.io_update_state(
           f_state=cur_state_file,f_eqdsk=cur_eqdsk_file,
           f_fastran='fastran.nc',time=timeStamp,relax=relax)

        #shutil.copyfile('fastran.nc',cur_fastran_file)

        #--- update plasma state files

        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            raise Exception, e

        #--- archive output files

        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            raise Exception, e

        return
    
    def finalize(self, timeStamp=0.0):

        return
    
