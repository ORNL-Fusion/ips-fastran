#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 esc component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *
from netCDF4 import *

from  component import Component

#-------------------
#--- zcode libraries
from Namelist import Namelist
import zefitutil,zefit
from zplasmastate import plasma_state_file

class esc(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- code entry

        print 'enter esc.step()'

        if (self.services == None) :
            print 'Error in esc.step () : No services'
            raise Exception('Error in esc.step (): No services')
        services = self.services

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- generate inefit

        zefit.io_input_from_state(
            cur_eqdsk_file,cur_state_file,cur_bc_file)

        #--- excutables

        try:
            esc_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            esc_bin = os.path.join(self.BIN_PATH, 'xesc')
        print esc_bin

        wgeqdsk_bin = os.path.join(self.BIN_PATH, 'wgeqdsk')
        print wgeqdsk_bin

        pstool_bin = os.path.join(self.BIN_PATH, 'pstool')
        print pstool_bin

        #--- run esc

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, esc_bin, logfile = 'xesc.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'esc'
           raise

        #--- run wgeqdsk

        logfile = open('wgeqdsk.log', 'w')
        retcode = subprocess.call([wgeqdsk_bin]
                      ,stdout=logfile,stderr=logfile,shell=True)
        logfile.close()

        if (retcode != 0):
           print 'Error executing ', 'wgeqdsk'
           raise
       
        #--- update local geqdsk state

        shutil.copyfile('geqdsk', cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        geq = zefitutil.readg("geqdsk") 
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']
        print 'r0 = ',r0
        print 'b0 = ',b0
        print 'ip = ',ip 

        shutil.copyfile(cur_state_file,"ps.nc")

        ps = plasma_state_file("ps.nc",r0=r0,b0=b0,ip=ip)
        j_tot = ps.dump_j_parallel()
        ps.close()

        logfile = open('pstool.log', 'w')
        retcode = subprocess.call([pstool_bin, "load", "geqdsk", "1.0d-6"],
                      stdout=logfile,stderr=logfile)
        logfile.close()
        if (retcode != 0):
           print 'Error executing ', pstool_bin
           raise

        ps = plasma_state_file("ps.nc",r0=r0,b0=b0,ip=ip)
        ps.load_j_parallel(j_tot)
        ps.close()

        shutil.copyfile("ps.nc",cur_state_file)

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
    
