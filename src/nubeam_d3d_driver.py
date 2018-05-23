#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 nubeam component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,os.path,shutil,pickle
import subprocess
from numpy import *

from  component import Component

# import zcode libraries
import Namelist
import znubeam

class nubeam_d3d_driver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        if (self.services == None) :
            print 'Error in nubeam.step() : No services'
            raise Exception('Error in nubeam.step(): No services')
        services = self.services

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e

        #--- get plasma state file names

        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- generate nubeam input

        instate = Namelist.Namelist(cur_instate_file)
        innubeam = Namelist.Namelist("innubeam")

        znubeam.io_write_input_files(instate,innubeam)
        shutil.copyfile(cur_eqdsk,"dnubeam_eqdsk.dat")

        #--- run nubeam init

        ncpu =  int(self.NPROC)
        nstep = innubeam["nubeam_run"]["nstep"][0]
        dt_nubeam = innubeam["nubeam_run"]["dt_nubeam"][0]
        print "ncpu  = ",ncpu
        print "dt    = ",dt_nubeam
        print "nstep = ",nstep

        print 'nubeam_init'
        workdir = services.get_working_dir()
        try:
            nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            nubeam_bin = os.path.join(self.BIN_PATH, 'dnubeam')
        print nubeam_bin

        task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, "init", logfile = 'dnubeam_init.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'nubeam init'
           raise

        for k in range(nstep):

            print 'nubeam_step %d'%k
            task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, "step", logfile = 'dnubeam_step.log')
            retcode = services.wait_task(task_id)

            if (retcode != 0):
               print 'Error executing ', 'nubeam step'
               raise
            shutil.copyfile("dnubeam_out.dat","dnubeam_out_%d.dat"%k)

        #--- update local instate

        znubeam.io_update_instate(
           nstep = nstep,
           navg  = 5,
           f_instate=cur_instate_file,
           f_outnubeam='dnubeam_out.dat')

        # shutil.copyfile("cur_state.cdf",cur_state_file)

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
    
