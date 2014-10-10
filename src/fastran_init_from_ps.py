#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran init component 
 JM
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
from zplasmastate import plasma_state_file,instate2ps

class fastran_init_from_ps (Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('fastran_init_from_ps.init() called')
        return

    def step(self, timeStamp):

        # --------------------------------------------------------------
        # entry

        print ('fastran_init_from_ps.step() called')

        if (self.services == None) :
            print 'Error in fastran_init_from_ps.step() : No services'
            raise Exception('Error in fastran_init_from_ps.step(): No services')
        services = self.services

        # --------------------------------------------------------------
        # check if this is a restart simulation

        try:
            mode = services.get_config_param('SIMULATION_MODE')
            print 'fastran_init: ',mode
        except Exception, e:
            print 'fastran_init: No SIMULATION_MODE variable in config file' \
                  ', NORMAL assumed', e
            mode='NORMAL'
            
        # --------------------------------------------------------------
        # RESTART simulation mode

        if mode == 'RESTART':
           #not implemented
           pass

        # --------------------------------------------------------------
        # NORMAL simulation mode
        
        elif mode == 'NORMAL':
        
            #----------------------------------------------------------
            #-- define run identifiers

            tokamak_id = services.get_config_param('TOKAMAK_ID')
            shot_number = services.get_config_param('SHOT_NUMBER')
            run_id = services.get_config_param('RUN_ID')

            print 'tokamak_id =', tokamak_id
            print 'shot_number =',shot_number
            print 'run_id =', run_id
    
            #----------------------------------------------------------
            #-- stage input files

            services.stage_input_files(self.INPUT_FILES)

            #----------------------------------------------------------
            #-- get plasma state file names

            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_bc_file = services.get_config_param('CURRENT_BC')

            #----------------------------------------------------------
            #-- get plasma state file names

            try:
                pstool_bin = os.path.join(self.BIN_PATH, self.BIN)
            except:
                pstool_bin = os.path.join(self.BIN_PATH, 'pstool')

            #----------------------------------------------------------
            #-- plasma state file

            try:
                shutil.copyfile("inps.nc",cur_state_file)
            except Exception,e:
                print e
                raise Exception('no plasma state file provided')

            #----------------------------------------------------------
            #-- plasma state file

            if "ingeqdsk" in self.INPUT_FILES:

                shutil.copyfile("ingeqdsk","geqdsk")

            else:

                print 'pstool dump geqdsk'

                shutil.copyfile(cur_state_file,"ps.nc")

                logfile = open('pstool_geqdsk.log', 'w')
                retcode = subprocess.call([pstool_bin, "dump", "geqdsk"],
                              stdout=logfile,stderr=logfile)
                if (retcode != 0):
                    print 'Error executing ', pstool_bin
                    raise
                logfile.close()

            #----------------------------------------------------------
            #-- boundary condition plasma state file

            geq = zefitutil.readg("geqdsk")

            inbc = Namelist()
            inbc["inbc"]["r0"] = [geq["rzero"]]
            inbc["inbc"]["b0"] = [geq["bcentr"]]

            inbc["inbc"]["ip"] = [geq["cpasma"]*1.0e-6]
            inbc["inbc"]["nbdry"] = [geq["nbdry"]]
            inbc["inbc"]["rbdry"] = geq["rbdry"]
            inbc["inbc"]["zbdry"] = geq["zbdry"]

            inbc.write(cur_bc_file) 

            #----------------------------------------------------------
            #-- copy to plasma state files

            shutil.copyfile("geqdsk",cur_eqdsk_file)
    
        # --------------------------------------------------------------
        # Update plasma state

        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            raise

        # --------------------------------------------------------------
        # archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def checkpoint(self, timeStamp=0.0):

        print 'fastran_init.checkpoint() called'
        
        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timeStamp, self.RESTART_FILES)
        
    def finalize(self, timeStamp=0.0):

        print 'fastran_init.finalize() called'

