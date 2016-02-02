#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran init component 
 -----------------------------------------------------------------------
"""

import sys,os,os.path,shutil,pickle
import subprocess
from numpy import *

from component import Component

from Namelist import Namelist
import zefit

import netCDF4
import zefitutil
import zplasmastate
import zinstate

class fastran_init (Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('fastran_init.init() called')
        return

    def step(self, timeStamp):

        #-- entry

        print ('fastran_init.step() called')
        services = self.services

        #-- check if this is a restart simulation

        try:
            mode = services.get_config_param('SIMULATION_MODE')
            print 'fastran_init: ',mode
        except Exception, e:
            print 'fastran_init: No SIMULATION_MODE variable in config file' \
                  ', NORMAL assumed', e
            mode = 'NORMAL'
            
        #-- RESTART simulation mode

        if mode == 'RESTART':
           #not implemented
           pass

        #-- NORMAL simulation mode
        
        elif mode == 'NORMAL':
        
            #-- define run identifiers

            tokamak_id = services.get_config_param('TOKAMAK_ID')
            shot_number = services.get_config_param('SHOT_NUMBER')
            run_id = services.get_config_param('RUN_ID')

            print 'tokamak_id =', tokamak_id
            print 'shot_number =',shot_number
            print 'run_id =', run_id
    
            #-- stage input files

            services.stage_input_files(self.INPUT_FILES)

            #-- get plasma state file names

            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_bc_file = services.get_config_param('CURRENT_BC')

            #-- initial plasma state file from instate

            ps = zplasmastate.zplasmastate('ips',1)

            if "instate" not in self.INPUT_FILES:
                raise Exception('no instate file provided')

            ps.init_from_instate (fn="instate",runid=run_id,shot=int(shot_number),t0=0.0,t1=0.0)

            #-- load intoric to plasma state

            if "intoric" in self.INPUT_FILES:
                print 'load intoric'
                #ps.load_intoric()

            #-- load innubeam to plasma state

            if "innubeam" in self.INPUT_FILES:

                print 'load innubeam'
                ps.load_innubeam()

            #-- initial equilibrium

            if "ingeqdsk" in self.INPUT_FILES:

                shutil.copyfile("ingeqdsk","geqdsk")

            else:

                print 'solve initial equilibrium'

                try:
                    gs_solver = self.INIT_EQ
                except:
                    gs_solver = 'esc'
                    pass

                if gs_solver == 'esc':

                    zefit.io_input_from_instate("instate")
                    esc_bin = os.path.join(self.BIN_PATH,'xesc')
                    wgeqdsk_bin = os.path.join(self.BIN_PATH, 'wgeqdsk')

                    logfile = open('xesc.log', 'w')
                    retcode = subprocess.call([esc_bin]
                                  ,stdout=logfile,stderr=logfile,shell=True)
                    if (retcode != 0):
                       print 'Error executing ', 'esc'
                       raise
                    logfile.close()

                    logfile = open('wgeqdsk.log', 'a')
                    retcode = subprocess.call([wgeqdsk_bin]
                                  ,stdout=logfile,stderr=logfile,shell=True)
                    if (retcode != 0):
                       print 'Error executing ', 'wgeqdsk'
                       raise
                    logfile.close()

                elif gs_solver == 'efit':

                    efit_bin = os.path.join(self.BIN_PATH, 'efitd90 129 129')

                    shot = 0
                    time = 0
                    kfile = "k%06d.%05d"%(shot,time)
                    args = "2\n 1\n "+kfile
                    command = 'echo \"%s\"'%args + ' | ' + efit_bin 

                    f=open("xefit","w")
                    f.write(command)
                    f.close()

                    zefit.io_input_init("instate")
                    zefit.fixbdry_kfile_init(shot,time,f_inefit="inefit")

                    cwd = services.get_working_dir()
                    task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
                    retcode = services.wait_task(task_id)
                    if (retcode != 0):
                       print 'Error executing ', 'efit'
                       raise
                    shutil.copyfile("g%06d.%05d"%(shot,time),"geqdsk") #<--------

            shutil.copyfile("geqdsk",cur_eqdsk_file)

            #-- load geqdsk to plasma state

            print 'pstool: load geqdsk'
            ps.load_geqdsk("geqdsk")

            #-- load instate to plasma state

            print 'instate2ps'
            zinstate.instate2ps("instate",ps)
            ps.store(cur_state_file)

            #-- boundary condition state file

            instate = Namelist("instate")["instate"]
            inbc = Namelist()
            inbc["inbc"]["r0"] = instate["r0"]
            inbc["inbc"]["b0"] = instate["b0"]
            inbc["inbc"]["ip"] = instate["ip"]
            inbc["inbc"]["nbdry"] = instate["nbdry"]
            inbc["inbc"]["rbdry"] = instate["rbdry"]
            inbc["inbc"]["zbdry"] = instate["zbdry"]
            inbc.write(cur_bc_file) 
    
        #-- Update plasma state

        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            raise

        #-- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def checkpoint(self, timeStamp=0.0):

        print 'fastran_init.checkpoint() called'
        
        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timeStamp, self.RESTART_FILES)
        
    def finalize(self, timeStamp=0.0):

        print 'fastran_init.finalize() called'

