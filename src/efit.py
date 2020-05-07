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
from inmetric_io import ps_to_inmetric
from efit_eqdsk import readg, writeg, scaleg

from numpy import *
import subprocess

class efit(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        #--- entry
        print('enter efit.init()')
        services = self.services

        #--- get shot and time
        self.ishot = int(services.get_config_param('SHOT_NUMBER'))
        self.itime = int(float(self.TIME_ID)) if hasattr(self, "TIME_ID") else int(services.get_config_param('TIME_ID'))

        print(20*'*')
        print(self.itime)

        #--- write efit run script
        scale_gs = int(getattr(self, "SCALE_GS", "0"))
        efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        kfile = "k%06d.%05d"%(self.ishot, self.itime)
        if scale_gs: kfile = kfile+"_s"
        args = "2\n 1\n "+kfile
        command = 'echo \"%s\"'%args + ' | ' + efit_bin
        f=open("xefit","w")
        f.write(command)
        f.close()

        #--- initial force free equilibrium
        init_run = int(getattr(self, "INIT_RUN", 0))
        print('init_run = ', init_run)

        if init_run:
            print ('INIT EFIT, force free')
            services.stage_plasma_state()
            self.initial_equilibrium()
            services.update_plasma_state()

    def step(self, timeid=0):
        #--- entry
        print ('enter efit.step()')
        services = self.services

        use_instate = int(getattr(self, "USE_INSTATE", "0"))

        #--- freeze
        ifreeze = int(getattr(self, "FREEZE", 10000))
        irefreeze = int(getattr(self, "REFREEZE", -10000))
        if timeid > ifreeze and timeid < irefreeze:
            print("EFIT FREEZE: timeid = %d, ifreeze = %d"%(timeid, ifreeze))
            return

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- get plasma state file names
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')

        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')
        elif ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        #--- generate inefit
        f_inefit = "inefit"
        mode = getattr(self, "PRESSURE", "kinetic")

        profile_relax = int(getattr(self, "PROFILE_RELAX", "0"))
        if timeid > 0 and profile_relax:
            shutil.copyfile(f_inefit,f_inefit+"0")

        betan_target = float(getattr(self, "BETAN_TARGET", "-1.0"))
        if ps_backend == 'PS' and use_instate == 0:
            efit_io.io_input_from_state(f_ps=cur_state_file, f_inbc=cur_bc_file, f_inefit=f_inefit, mode=mode, betan_target=betan_target)
            nrho = 257
        elif ps_backend == 'INSTATE' or use_instate == 1:
            efit_io.io_input_from_instate(f_instate=cur_instate_file, f_inefit=f_inefit, mode=mode)
            instate = Namelist(cur_instate_file,"r")
            nrho = instate["instate"]["nrho"][0]

        if timeid > 0 and profile_relax:
            inefit0 = Namelist(f_inefit+"0")
            inefit = Namelist(f_inefit)
            for key in ["press", "jpar"]:
                inefit['inefit'][key] = 0.5*array(inefit0['inefit'][key])+0.5*array(inefit['inefit'][key])
            inefit.write(f_inefit)

        #--- scaled GS
        scale_gs = int(getattr(self, "SCALE_GS", "0"))
        R0_scale = float(getattr(self, "R0_scale", "0"))
        B0_scale = float(getattr(self, "B0_scale", "0"))

        if scale_gs:
            if ps_backend == 'PS':
                inbc = Namelist(cur_bc_file, 'r')
                Rs = inbc["inbc"]["r0"][0]/R0_scale
                Bs = inbc["inbc"]["b0"][0]/B0_scale
            else:
                Rs = instate["instate"]["r0"][0]/R0_scale
                Bs = instate["instate"]["b0"][0]/B0_scale
        else:
            Rs = 1.0
            Bs = 1.0

        #--- topology
        topology = getattr(self, "TOPOLOGY", "")
        print('topology =', topology)

        error = float(getattr(self, "ERROR", "1.0e-4"))
        print('error =', error)

        #--- clean up
        try:
           for f in glob.glob("g000000.?????"): os.remove(f)
           print ('g00000.* removed')
        except:
           pass

        #--- run efit, force-free
        init_run = int(getattr(self, "INIT_RUN_STEP", 0))
        print ('init_run = ', init_run)

        cwd = services.get_working_dir()

        niter = int(getattr(self, "NITER", "5"))

        for k in range(niter):
            print ("*********************")
            print ("*** generate kfile %d"%k)

            t0 = timer.time()

            if k == 0 and init_run:
                efit_io.fixbdry_kfile_init(self.ishot, self.itime, f_inefit=f_inefit)
                if scale_gs: efit_io.scale_kfile(self.ishot, self.itime, Rs=Rs, Bs=Bs)
                iconv = 0

            else:
                if k > 0:
                    relax = 1
                    print ('apply underrelaxation')
                else:
                    relax = 0

                ps = plasmastate('ips',1)

                geqdsk = readg(cur_eqdsk_file)

                r0  = geqdsk["rzero" ]
                b0  = abs(geqdsk["bcentr"])
                ip  = geqdsk['cpasma']
                print ('r0 = ',r0)
                print ('b0 = ',b0)
                print ('ip = ',ip)

                ps.init_from_geqdsk (cur_eqdsk_file,nrho=nrho,nth=101)
                inmetric = ps_to_inmetric(ps, r0, b0, ip)

                iconv = efit_io.fixbdry_kfile(self.ishot, self.itime, Namelist(f_inefit,"r"), inmetric["inmetric"], relax=relax, topology=topology, error=error)
                if scale_gs: efit_io.scale_kfile(self.ishot, self.itime, Rs=Rs, Bs=Bs)

            if iconv:
                print ('converged')
                break

            t1 = timer.time()

            print ("run efit")

            if int(getattr(self,'SERIAL','0')) == 1:
                print ('efit, subprocess')
                logfile = open('efit.log', 'w')
                retcode = subprocess.call(["sh xefit"]
                              ,stdout=logfile,stderr=logfile,shell=True)
                logfile.close()
            else:
                task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
                retcode = services.wait_task(task_id)

            if (retcode != 0):
               raise Exception('Error executing efit')

            t2 = timer.time()

            print ('**** %6.3f %6.3f'%(t1-t0, t2-t1))

            if scale_gs:
                shutil.copyfile("g%06d.%05d"%(self.ishot, self.itime), "g%06d.%05d_s"%(self.ishot, self.itime))
                g = readg("g%06d.%05d"%(self.ishot, self.itime))
                scaleg(g, R0=Rs, B0=Bs)
                writeg(g,self.ishot, self.itime, 129)

            #--- update local geqdsk state
            if cur_eqdsk_file != "g%06d.%05d"%(self.ishot, self.itime):
                shutil.copyfile("g%06d.%05d"%(self.ishot, self.itime), cur_eqdsk_file)

        #--- load geqdsk to plasma state file
        if ps_backend == 'PS':
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)

            ps_geqdsk = plasmastate('ips',1)
            ps_geqdsk.init_from_geqdsk (cur_eqdsk_file, nrho=ps["nrho"], nth=101)
            ps["vol"][:] = ps_geqdsk["vol"][:]
            ps["g_eq"][:] = ps_geqdsk["g_eq"][:]

            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':
            instate["inmetric"] = inmetric["inmetric"]
            instate.write(cur_instate_file)

        #--- update plasma state files
        services.update_plasma_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid):
        print ('enter efit.step()')
        services = self.services

        dir_state = services.get_config_param('PLASMA_STATE_WORK_DIR')
        for header in ['k','a']:
            f = "%s%06d.%05d"%(header,self.ishot,self.itime)
            shutil.copyfile(f, os.path.join(dir_state,f))

    def initial_equilibrium(self):
        services = self.services

        ishot = int(services.get_config_param('SHOT_NUMBER'))
        itime = int(services.get_config_param('TIME_ID'))

        scale_gs = int(getattr(self, "SCALE_GS", "0"))
        R0_scale = float(getattr(self, 'R0_scale', '0'))
        B0_scale = float(getattr(self, 'B0_scale', '0'))

        #--- initial force free equilibrium

        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')
            print (cur_bc_file)
            print (Namelist(cur_bc_file))
            instate = Namelist(cur_bc_file)["inbc"]
        elif ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')
            instate = Namelist(cur_instate_file)["instate"]

        nrho = 101
        inefit = Namelist()
        inefit["inefit"]["ip"   ] = [instate["ip"][0]*1.0e6]
        inefit["inefit"]["r0"   ] = instate["r0"]
        inefit["inefit"]["b0"   ] = instate["b0"]
        inefit["inefit"]["nrho" ] = [nrho]
        inefit["inefit"]["rho"  ] = linspace(0,1.0,nrho)
        inefit["inefit"]["press"] = nrho*[0.0]
        inefit["inefit"]["jpar" ] = nrho*[0.0]
        inefit["inefit"]["nlim" ] = instate["nlim" ]
        inefit["inefit"]["rlim" ] = instate["rlim" ]
        inefit["inefit"]["zlim" ] = instate["zlim" ]
        inefit["inefit"]["nbdry"] = instate["nbdry"]
        inefit["inefit"]["rbdry"] = instate["rbdry"]
        inefit["inefit"]["zbdry"] = instate["zbdry"]
        inefit.write("inefit")

        efit_io.fixbdry_kfile_init(ishot, itime, f_inefit="inefit")

        if scale_gs:
            print('scaled ',R0_scale, B0_scale)
            Rs = instate["r0"][0]/R0_scale
            Bs = instate["b0"][0]/B0_scale
            efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)

        f=open("log.efit", "w")
        subprocess.call(["sh", "xefit"], stdout=f)

        if scale_gs:
            g = readg("g%06d.%05d"%(ishot, itime))
            scaleg(g, R0=Rs, B0=Bs)
            writeg(g, ishot, itime, 129)

        if cur_eqdsk_file != "g%06d.%05d"%(ishot, itime):
            shutil.copyfile("g%06d.%05d"%(self.ishot, self.itime), cur_eqdsk_file)

        if ps_backend == 'PS':
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps_geqdsk = plasmastate('ips',1)
            ps_geqdsk.init_from_geqdsk (cur_eqdsk_file, nrho=ps["nrho"], nth=101)
            ps["vol"][:] = ps_geqdsk["vol"][:]
            ps["g_eq"][:] = ps_geqdsk["g_eq"][:]
            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':
            instate["inmetric"] = inmetric["inmetric"]
            instate.write(cur_instate_file)
