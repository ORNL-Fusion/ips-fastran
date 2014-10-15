#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 efit component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from  component import Component

#-------------------
#--- zcode libraries
import zefit, zefitutil
from zplasmastate import plasma_state_file

class efit(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        print 'enter efit.step()'

        if (self.services == None) :
            print 'Error in efit.step () : No services'
            raise Exception('Error in efit.step (): No services')
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

        #--- run efit

        try:
            efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            efit_bin = os.path.join(self.BIN, 'efitd90 129 129')

        niter = 5
        shot=0
        time=0
        f_inefit="inefit"

        logfile = open('efit.log', 'w')
        logfile.close()
        for k in range(niter):

            print "generate kfile %d"%k
            if k == 0:
                zefit.fixbdry_kfile_init(shot,time,f_inefit)
            else:
                zefit.fixbdry_kfile(shot,time,f_inefit)
    
            print "run efit"
            kfile = "k000000.00000"
            args = "2\n 1\n "+kfile
            command = 'echo \"%s\"'%args + ' | ' + efit_bin 
            logfile = open('efit.log', 'a')
            retcode = subprocess.call([command]
                          ,stdout=logfile,stderr=logfile,shell=True)
            if (retcode != 0):
               print 'Error executing ', 'efit'
               raise
            logfile.close()

        #--- update local geqdsk state

        shutil.copyfile("g000000.00000", cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        shutil.copyfile("g000000.00000", "geqdsk")
        shutil.copyfile(cur_state_file,"ps.nc")

        try:
            FASTRAN_ROOT = os.environ["FASTRAN_ROOT"]
            DIR_BIN = os.path.join(FASTRAN_ROOT,"bin")
            pstool_bin = os.path.join(DIR_BIN,"pstool")
        except:
            pstool_bin = 'pstool'

        geq = zefitutil.readg("geqdsk") 
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']
        print 'r0 = ',r0
        print 'b0 = ',b0
        print 'ip = ',ip 

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
    
