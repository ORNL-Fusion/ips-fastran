#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 efit component 
 -----------------------------------------------------------------------
"""

import sys,os,shutil,glob
import subprocess
from numpy import *
import time as timer

#--- ips framework
from  component import Component

#--- zcode libraries
import zefit, zefitutil
import zplasmastate

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
        try:
            time = int(services.get_config_param('TIME_ID'))
        except:
            time = 0

        #--- generate inefit

        try:
            mode = self.PRESSURE
        except:
            mode = 'kinetic'
        print 'mode = ', mode

        zefit.io_input_from_state(cur_state_file,cur_bc_file,mode=mode)

        #--- clean up

        try: 
           files = glob.glob("g000000.?????")
           for file in files: os.remove(file)
           print 'removed'
        except:
           pass

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

            print "*** generate kfile %d"%k

            t0 = timer.time()
 
            if k == 0:
                iconv = zefit.fixbdry_kfile_init(shot,time,f_inefit)
            else:
                if k > kinit: 
                   relax = 1
                   print 'apply underrelaxation'
                else: relax = 0
                iconv = zefit.fixbdry_kfile(shot,time,f_inefit,relax=relax)
         
            if iconv:
               print 'converged'
               break

            t1 = timer.time()

            print "run efit"

            task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
            retcode = services.wait_task(task_id)
            if (retcode != 0):
               print 'Error executing ', 'efit'
               raise
            shutil.copyfile("g%06d.%05d"%(shot,time),"g000000.%05d"%(k)) #<--------

            t2 = timer.time()

            print '**** %6.3f %6.3f'%(t1-t0,t2-t1)

        #--- update local geqdsk state

        shutil.copyfile("g%06d.%05d"%(shot,time), cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        geq = zefitutil.readg(cur_eqdsk_file) 
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']
        print 'r0 = ',r0
        print 'b0 = ',b0
        print 'ip = ',ip 

        ps = zplasmastate.zplasmastate('ips',1)
        ps.read(cur_state_file)

        j_tot = 1.e-6*ps.dump_j_parallel(ps["rho"],"rho_eq","curt",r0,b0,tot=True)
        ps.load_geqdsk(cur_eqdsk_file)
        ps.load_j_parallel(ps["rho"],j_tot,"rho_eq","curt",r0,b0,tot=True)

        ps.store(cur_state_file)

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
    
