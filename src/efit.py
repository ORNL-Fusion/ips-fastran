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

    def init(self, timeStamp=0):

        return

    def step(self, timeStamp=0):

        #--- entry

        print 'enter efit.step()'

        services = self.services

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- get shot and time

        try:
            shot = int(services.get_config_param('SHOT_NUMBER'))
        except:
            shot = 0
        time = 0

        #--- generate inefit

        zefit.io_input_from_state(cur_state_file,cur_bc_file)

        #--- run efit

        try:
            efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            efit_bin = os.path.join(self.BIN, 'efitd90 129 129')

        try:
            niter = int(self.NITER)
        except:
            niter = 5
        try:
            kinit = int(self.RESTART)
        except:
            kinit = 0

        shutil.copyfile(cur_eqdsk_file,"g%06d.%05d"%(shot,time))

        f_inefit="inefit"

        kfile = "k%06d.%05d"%(shot,time)
        args = "2\n 1\n "+kfile
        command = 'echo \"%s\"'%args + ' | ' + efit_bin 

        f=open("xefit","w")
        f.write(command)
        f.close()

        cwd = services.get_working_dir()

        for k in range(kinit,kinit+niter):

            print "generate kfile %d"%k
            if k == 0:
                zefit.fixbdry_kfile_init(shot,time,f_inefit)
            else:
                zefit.fixbdry_kfile(shot,time,f_inefit)
    
            print "run efit"

            task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
            retcode = services.wait_task(task_id)
            if (retcode != 0):
               print 'Error executing ', 'efit'
               raise

        #--- update local geqdsk state

        shutil.copyfile("g%06d.%05d"%(shot,time), cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        shutil.copyfile("g%06d.%05d"%(shot,time), "geqdsk")
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
    
    def finalize(self, timeStamp=0):
        return
    
