#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 efit component
 -----------------------------------------------------------------------
"""

import os
import shutil
import glob
import time as timer

from component import Component

from Namelist import Namelist
import efit_io
from plasmastate import plasmastate
from zinmetric import zinmetric
from efit_eqdsk import readg, writeg,scaleg

from numpy import *

class efit(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        #--- entry

        print 'enter efit.init()'
        services = self.services

        #--- get shot and time

        self.ishot = int(services.get_config_param('SHOT_NUMBER'))
        self.itime = int(float(self.TIME_ID)) if hasattr(self,"TIME_ID") else int(services.get_config_param('TIME_ID'))

        print 20*'*'
        print self.itime

        #--- for scaled GS

        R0_scale = float(self.R0_scale) if hasattr(self,"R0_scale") else 0
        B0_scale = float(self.B0_scale) if hasattr(self,"B0_scale") else 0

        if R0_scale > 0 and B0_scale > 0:
            scaled_gs = 1
        else:
            scaled_gs = 0

        #--- write efit run script

        efit_bin = os.path.join(self.BIN_PATH, self.BIN)

        if scaled_gs:
            kfile = "k%06d.%05d_s"%(self.ishot,self.itime)
        else:
            kfile = "k%06d.%05d"%(self.ishot,self.itime)

        args = "2\n 1\n "+kfile
        command = 'echo \"%s\"'%args + ' | ' + efit_bin

        f=open("xefit","w")
        f.write(command)
        f.close()

        #--- initial force free equilibrium

        ps_backend = getattr(self,'PS_BACKEND','PS')

        init_run = int(getattr(self,"INIT_RUN",0))
        print 'init_run = ',init_run

        if init_run:

            print 'INIT EFIT, force free'

            services.stage_plasma_state()

            cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')

            services.stage_input_files(self.INPUT_FILES)
            if ps_backend == 'PS':
                efit_io.io_input_from_state(cur_state_file,cur_bc_file,f_inefit="inefit",mode='total')
            elif ps_backend == 'INSTATE':
                #efit_io.io_input_init(cur_bc_file)
                cur_instate_file = services.get_config_param('CURRENT_INSTATE')
                efit_io.io_input_init(cur_instate_file)

            efit_io.fixbdry_kfile_init(self.ishot,self.itime,f_inefit="inefit")

            if scaled_gs:
                inbc = Namelist(cur_bc_file,'r')
                Rs = inbc["inbc"]["r0"][0]/R0_scale
                Bs = inbc["inbc"]["b0"][0]/B0_scale
                efit_io.scale_kfile(self.ishot,self.itime,Rs=Rs,Bs=Bs)

            cwd = services.get_working_dir()
            task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
            retcode = services.wait_task(task_id)

            if (retcode != 0):
                raise Exception("Error executing efit")

            if scaled_gs:
                g = readg("g%06d.%05d"%(self.ishot,self.itime))
                scaleg(g,R0=Rs,B0=Bs)
                writeg(g,self.ishot,self.itime,129)

            if cur_eqdsk_file != "g%06d.%05d"%(self.ishot,self.itime):
                shutil.copyfile("g%06d.%05d"%(self.ishot,self.itime), cur_eqdsk_file)

            if ps_backend == 'PS':
                print '************* LOAD'
                ps = plasmastate('ips',1)
                ps.read(cur_state_file)
                ps.load_geqdsk(cur_eqdsk_file)
                ps.store(cur_state_file)
                services.update_plasma_state()
                print '************* END'

        return

    def step(self, timeStamp=0):

        #--- entry

        print 'enter efit.step()'

        services = self.services

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file names

        ps_backend = getattr(self,'PS_BACKEND','PS')

        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')
        elif ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        #--- generate inefit

        f_inefit = "inefit"
        mode = self.PRESSURE if hasattr(self,"PRESSURE") else "kinetic"

        profile_relax = int(getattr(self,"PROFILE_RELAX","1"))

        if timeStamp > 0 and profile_relax:
            shutil.copyfile(f_inefit,f_inefit+"0")

        betan_target = float(getattr(self,"BETAN_TARGET","-1.0"))
        if ps_backend == 'PS':
            efit_io.io_input_from_state(f_ps=cur_state_file,f_inbc=cur_bc_file,f_inefit=f_inefit,mode=mode,betan_target=betan_target)
            nrho = 257
        elif ps_backend == 'INSTATE':
            efit_io.io_input_from_instate(f_instate=cur_instate_file, f_inefit=f_inefit, mode=mode)
            instate = Namelist(cur_instate_file,"r")
            nrho = instate["instate"]["nrho"][0]

        if timeStamp > 0 and profile_relax:
            inefit0 = Namelist(f_inefit+"0")
            inefit = Namelist(f_inefit)
            for key in ["press","jpar"]:
                inefit['inefit'][key] = 0.5*array(inefit0['inefit'][key])+0.5*array(inefit['inefit'][key])
            inefit.write(f_inefit)

        #--- scaled GS

        R0_scale = float(self.R0_scale) if hasattr(self,"R0_scale") else 0
        B0_scale = float(self.B0_scale) if hasattr(self,"B0_scale") else 0

        if R0_scale > 0 and B0_scale > 0:
            scaled_gs = 1
        else:
            scaled_gs = 0

        if scaled_gs:
            if ps_backend == 'PS':
                inbc = Namelist(cur_bc_file,'r')
                Rs = inbc["inbc"]["r0"][0]/R0_scale
                Bs = inbc["inbc"]["b0"][0]/B0_scale
            else:
                Rs = instate["instate"]["r0"][0]/R0_scale
                Bs = instate["instate"]["b0"][0]/B0_scale

        else:
            Rs = 1.0
            Bs = 1.0

        topology = getattr(self,"TOPOLOGY","")
        print 'topology =', topology

        error = float(getattr(self,"ERROR","1.0e-4"))
        print 'error =', error

        #--- clean up

        try:
           for f in glob.glob("g000000.?????"): os.remove(f)
           print 'g00000.* removed'
        except:
           pass

        #--- run efit, force-free

        init_run = int(getattr(self,"INIT_RUN_STEP",0))
        print 'init_run = ',init_run

        cwd = services.get_working_dir()

        niter = int(self.NITER) if hasattr(self,"NITER") else 5

        for k in range(niter):

            print "*********************"
            print "*** generate kfile %d"%k

            t0 = timer.time()

            #if k == 0 and timeStamp ==-100:
            if k == 0 and init_run:

                efit_io.fixbdry_kfile_init(self.ishot,self.itime,f_inefit=f_inefit)
                if scaled_gs: efit_io.scale_kfile(self.ishot,self.itime,Rs=Rs,Bs=Bs)
                iconv = 0

            else:

                if k > 0:
                    relax = 1
                    print 'apply underrelaxation'
                else:
                    relax = 0

                ps = plasmastate('ips',1)

                geqdsk = readg(cur_eqdsk_file)

                r0  = geqdsk["rzero" ]
                b0  = abs(geqdsk["bcentr"])
                ip  = geqdsk['cpasma']
                print 'r0 = ',r0
                print 'b0 = ',b0
                print 'ip = ',ip

                ps.init_from_geqdsk (cur_eqdsk_file,nrho=nrho,nth=101)
                inmetric = zinmetric(ps,r0,b0,ip)

                iconv = efit_io.fixbdry_kfile(self.ishot,self.itime,Namelist(f_inefit,"r"),inmetric["inmetric"],relax=relax,topology=topology,error=error)
                if scaled_gs: efit_io.scale_kfile(self.ishot,self.itime,Rs=Rs,Bs=Bs)

            if iconv:
                print 'converged'
                break

            t1 = timer.time()

            print "run efit"

            task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
            retcode = services.wait_task(task_id)
            if (retcode != 0):
               raise Exception('Error executing efit')
           #shutil.copyfile("g%06d.%05d"%(self.ishot,self.itime),"g000000.%05d"%(k))

            t2 = timer.time()

            print '**** %6.3f %6.3f'%(t1-t0,t2-t1)

            if scaled_gs:
                g = readg("g%06d.%05d"%(self.ishot,self.itime))
                scaleg(g,R0=Rs,B0=Bs)
                writeg(g,self.ishot,self.itime,129)

            #--- update local geqdsk state

            if cur_eqdsk_file != "g%06d.%05d"%(self.ishot,self.itime):
                shutil.copyfile("g%06d.%05d"%(self.ishot,self.itime), cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        if ps_backend == 'PS':

            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':

            instate["inmetric"] = inmetric["inmetric"]
            instate.write(cur_instate_file)

        #--- update plasma state files

        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def finalize(self, timeStamp=0):

        services = self.services

        dir_state = services.get_config_param('PLASMA_STATE_WORK_DIR')
        for header in ['k','a']:
            f = "%s%06d.%05d"%(header,self.ishot,self.itime)
            shutil.copyfile(f, os.path.join(dir_state,f))
        print 'HAHA'

        return
