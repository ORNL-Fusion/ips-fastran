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
from inspect import currentframe, getframeinfo

class fastran_init (Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('fastran_init.init() called')
        return

    def step(self, timeStamp):

        # --------------------------------------------------------------
        # entry

        print ('fastran_init.step() called')

        if (self.services == None) :
            print 'Error in fastran_init.step() : No services'
            raise Exception('Error in fastran_init.step(): No services')
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
            #-- get pstool excutable name

            try:
                pstool_bin = os.path.join(self.BIN_PATH, self.BIN)
            except:
                pstool_bin = os.path.join(self.BIN_PATH, 'pstool')
            print 'pstool_bin path:'
            print pstool_bin


            #----------------------------------------------------------
            #-- initial plasma state file from instate

            if "instate" not in self.INPUT_FILES:
                raise Exception('no instate file provided')

            instate = Namelist("instate")["instate"]
            nrho = instate["nrho"][0]
            print "nrho = ", nrho

            inps = Namelist()

            inps["inps"]["global_label"] = ['fastran_scenario']
            inps["inps"]["tokamak_id"  ] = [tokamak_id]
            inps["inps"]["runid"       ] = [run_id]

            inps["inps"]["shot_number" ] = [int(shot_number)]
            inps["inps"]["time"        ] = [0.0]

            inps["inps"]["nspec_ion"   ] = instate["n_ion"] 
            inps["inps"]["nspec_imp"   ] = instate["n_imp"]
            inps["inps"]["nspec_beam"  ] = [1]
            inps["inps"]["nspec_fusion"] = [1]
            inps["inps"]["nspec_rfmin" ] = [1]
            inps["inps"]["nspec_gas"   ] = instate["n_ion"] 
         
            inps["inps"]["z_ion"       ] = instate["z_ion"]
            inps["inps"]["a_ion"       ] = instate["a_ion"]
            inps["inps"]["z_imp"       ] = instate["z_imp"]
            inps["inps"]["a_imp"       ] = instate["a_imp"]
            inps["inps"]["z_beam"      ] = [1]
            inps["inps"]["a_beam"      ] = [2]
            inps["inps"]["z_fusion"    ] = [2]
            inps["inps"]["a_fusion"    ] = [4]
            inps["inps"]["z_rfmin"     ] = [2]
            inps["inps"]["a_rfmin"     ] = [3]
            inps["inps"]["z_gas"       ] = instate["z_ion"]
            inps["inps"]["a_gas"       ] = instate["a_ion"]

            inps["inps"]["nrho"        ] = [nrho]
            inps["inps"]["nrho_eq"     ] = [nrho]
            inps["inps"]["nth_eq"      ] = [101]
            inps["inps"]["nrho_eq_geo" ] = [nrho]
            inps["inps"]["nrho_gas"    ] = [nrho]
            inps["inps"]["nrho_nbi"    ] = [nrho]
            inps["inps"]["nrho_ecrf"   ] = [nrho]
            inps["inps"]["nrho_icrf"   ] = [nrho]
            inps["inps"]["nrho_fus"    ] = [nrho]
            inps["inps"]["nrho_anom"   ] = [nrho]

            inps.write("inps")

            print 'pstool init'

            logfile=open("pstool_init.log","w")
            retcode = subprocess.call([pstool_bin, "init"],
                          stdout=logfile,stderr=logfile)
            logfile.close()
            if (retcode != 0):
               raise Exception('Error executing ', pstool_bin)

            #----------------------------------------------------------
            #-- load intoric to plasma state

            if "intoric" in self.INPUT_FILES:

                print 'pstool load inrf'
    
                logfile=open("pstool_intoric.log","w")
                retcode = subprocess.call([pstool_bin, "load","intoric"],
                              stdout=logfile,stderr=logfile)
                logfile.close()
                if (retcode != 0):
                   raise Exception('Error executing ', pstool_bin)

            #----------------------------------------------------------
            #-- load innubeam to plasma state

            if "innubeam" in self.INPUT_FILES:
                print 'pstool load innubeam'
                logfile=open("pstool_innubeam.log","w")
                retcode = subprocess.call([pstool_bin, "load", "innubeam"],
                              stdout=logfile,stderr=logfile)
                logfile.close()
                if (retcode != 0):
                   print 'Error executing ', pstool_bin
                   raise

            #----------------------------------------------------------
            #-- initial equilibrium

            if "ingeqdsk" in self.INPUT_FILES:

                shutil.copyfile("ingeqdsk","geqdsk")

            else:

                print 'solve initial equilibrium'

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

            #----------------------------------------------------------
            #-- load geqdsk to plasma state

            print 'pstool load geqdsk'

            logfile = open('pstool_geqdsk.log', 'w')
            retcode = subprocess.call([pstool_bin, "load", "geqdsk", "1.0d-6"],
                          stdout=logfile,stderr=logfile)
            if (retcode != 0):
               print 'Error executing ', pstool_bin
               raise
            logfile.close()

            #----------------------------------------------------------
            #-- load instate to plasma state

            print 'instate2ps'

            geq = zefitutil.readg("geqdsk") 
            r0  = geq["rzero" ]
            b0  = abs(geq["bcentr"])
            ip  = geq['cpasma']

            ps = plasma_state_file("ps.nc",r0=r0,b0=b0,ip=ip)
            instate2ps(instate,ps)
            ps.close()

            #----------------------------------------------------------
            #-- boundary condition plasma state file

            inbc = Namelist()
            inbc["inbc"]["r0"] = instate["r0"]
            inbc["inbc"]["b0"] = instate["b0"]
            inbc["inbc"]["ip"] = instate["ip"]
            inbc["inbc"]["nbdry"] = instate["nbdry"]
            inbc["inbc"]["rbdry"] = instate["rbdry"]
            inbc["inbc"]["zbdry"] = instate["zbdry"]
            inbc.write(cur_bc_file) 
    
            #----------------------------------------------------------
            #-- copy to plasma state files

            shutil.copyfile("ps.nc",cur_state_file)
            shutil.copyfile("geqdsk",cur_eqdsk_file)

            #----------------------------------------------------------
            #-- copy current plasma state to prior state and next state
            if "PRIOR_STATE" in self.PLASMA_STATE_FILES:
                try:
                    prior_state_file = services.get_config_param('PRIOR_STATE')
                    shutil.copyfile(cur_state_file, prior_state_file)
                except Exception, e:
                    frameinfo = getframeinfo(currentframe())
                    print framemeinfo.filename, frameinfo.lineno
                    print 'No PRIOR_STATE file ', e
                    raise

            if "NEXT_STATE" in self.PLASMA_STATE_FILES:
                try:
                    next_state_file = services.get_config_param('NEXT_STATE')
                    shutil.copyfile(cur_state_file, next_state_file)
                except Exception, e:
                    frameinfo = getframeinfo(currentframe())
                    print framemeinfo.filename, frameinfo.lineno
                    print 'No NEXT_STATE file ', e
                    raise

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

