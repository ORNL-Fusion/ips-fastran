"""
 -----------------------------------------------------------------------
 efit component
 -----------------------------------------------------------------------
"""
import os
import re
import shutil
import glob
import subprocess
import time as timer
from numpy import *
from Namelist import Namelist
from fastran.equilibrium import efit_io
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver.inmetric_io import ps_to_inmetric
from fastran.equilibrium.efit_eqdsk import readg, writeg, scaleg
from ipsframework import Component

class efit(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def _get_shot_time(self):
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(float(self.TIME_ID)) if hasattr(self, "TIME_ID") else int(self.services.get_config_param('TIME_ID'))
        return ishot, itime 

    def _write_xefit(self, ishot, itime):
        scale_gs = int(getattr(self, "SCALE_GS", "0"))
        efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        kfile = "k%06d.%05d"%(ishot, itime)
        if scale_gs: kfile = kfile + "_s"
        args = "2\n 1\n " + kfile
        command = 'echo \"%s\"'%args + ' | ' + efit_bin
        f=open("xefit", "w")
        f.write(command)
        f.close()

    def init(self, timeid=0):
        print('>>> efit.init() started')
        services = self.services

        #--- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid

        print('efit time = {}'.format(itime))

        #--- write efit run script
        self._write_xefit(ishot, itime)

        #--- initial force free equilibrium
        init_run = int(getattr(self, "INIT_RUN", 0))
        print('init_run = ', init_run)

        if init_run:
            print('force free')
            services.stage_state()
            self.initial_equilibrium(ishot, itime)
            services.update_state()

    def step(self, timeid=0):
        #--- entry
        print ('>>> efit.step() started')
        services = self.services

        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid

        #--- freeze
        ifreeze = int(getattr(self, "FREEZE", 10000))
        iresume = int(getattr(self, "RESUME", -10000))
        if timeid > ifreeze and timeid < iresume:
            print("EFIT FREEZE: timeid = {}, ifreeze = {}, iresume = {}"%(timeid, ifreeze, iresume))
            return

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')

        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')
        elif ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        use_instate = int(getattr(self, "USE_INSTATE", "0"))

        #--- generate inefit
        f_inefit = "inefit"
        mode = getattr(self, "PRESSURE", "kinetic")
        print ('mode =', mode)

        profile_relax = int(getattr(self, "PROFILE_RELAX", "0"))
        if timeid >= 0 and profile_relax:
            shutil.copyfile(f_inefit, f_inefit+"0")

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
                inefit['inefit'][key] = 0.5*array(inefit0['inefit'][key]) + 0.5*array(inefit['inefit'][key])
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

        max_error_efit = float(getattr(self, "MAX_EFIT_ERROR", "0.1"))
        print('max_error_efit =', max_error_efit)

        #--- clean up
        try:
           for f in glob.glob("g000000.?????"): os.remove(f)
           print ('g00000.* removed')
        except:
           pass

        #--- write efit run script
        self._write_xefit(ishot, itime)

        #--- run efit, force-free
        init_run_step = int(getattr(self, "INIT_RUN_STEP", 0))
        print ('init_run_step =', init_run_step)

        cwd = services.get_working_dir()

        niter = int(getattr(self, "NITER", "5"))
        print ('niter = ', niter)

        pat = re.compile("\s*it=")
        for k in range(niter):
            print ("*********************")
            print ("*** generate kfile %d"%k)

            t0 = timer.time()

            if k == 0 and init_run_step:
                efit_io.fixbdry_kfile_init(ishot, itime, f_inefit=f_inefit)
                if scale_gs: efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)
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

                try:
                    ps.init_from_geqdsk (cur_eqdsk_file, nrho=nrho, nth=101)
                except:
                    raise Exception("erro in init_from_geqdsk")
                if k > 0:
                    inmetric_prev = inmetric
                inmetric = ps_to_inmetric(ps, r0, b0, ip)
                inmetric_keys = [
                  'volp', 'ipol', 'g11', 'g22', 'g33', 'gradrho', 'area', 'rminor', 'rmajor',
                  'shift', 'kappa', 'delta', 'pmhd', 'qmhd',
                  'er', 'nc1', 'hfac1', 'hfac2', 'psi', "vol", "gr2i", "bp2"]
                if k > 0:
                    print ("inmtric under-relax")
                    for key in inmetric_keys:
                        inmetric["inmetric"][key] = 0.5*array(inmetric_prev["inmetric"][key]) + 0.5*array(inmetric["inmetric"][key])

                iconv = efit_io.fixbdry_kfile(ishot, itime, Namelist(f_inefit, "r"), inmetric["inmetric"], relax=relax, topology=topology, error=error)
                if scale_gs: efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)

            if iconv:
                print ('converged')
                break

            t1 = timer.time()

            print("run efit")

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

            print('elapsed time = %6.3f %6.3f'%(t1-t0, t2-t1))

            error_efit = []
            logfile = open("efit.log", "r")
            lines = logfile.readlines()
            for line in lines:
                if pat.search(line):
                    try:
                        e = float(line.split("err=")[1].split()[0])
                    except:
                        e = 0.
                    error_efit.append(e)
            logfile.close()
            print('error =', error_efit[-1])
            if abs(error_efit[-1]) > max_error_efit:
               raise Exception("EFIT not converged")

            if scale_gs:
                shutil.copyfile("g%06d.%05d"%(ishot, itime), "g%06d.%05d_s"%(ishot, itime))
                g = readg("g%06d.%05d"%(ishot, itime))
                scaleg(g, R0=Rs, B0=Bs)
                writeg(g,ishot, itime, 129)

            #--- update local geqdsk state
            if cur_eqdsk_file != "g%06d.%05d"%(ishot, itime):
                shutil.copyfile("g%06d.%05d"%(ishot, itime), cur_eqdsk_file)

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
        services.update_state()

        #--- archive output files
        if self.OUTPUT_FILES:
            services.stage_output_files(timeid, self.OUTPUT_FILES)
        else:
            output_files = ' '.join(["{0}{1:06d}.{2:05d}".format(header, ishot, itime) for header in ["g", "a", "k", "m"]])
            output_files += " inefit"
            print("stage_output_files: ", output_files)
            services.stage_output_files(timeid,  output_files)

    def finalize(self, timeid):
        print ('>>> efit.finalize() started')
        services = self.services

        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid

        dir_state = services.get_config_param('STATE_WORK_DIR')
        for header in ['k', 'a']:
            f = "%s%06d.%05d"%(header, ishot, itime)
            shutil.copyfile(f, os.path.join(dir_state, f))

    def initial_equilibrium(self, ishot, itime):
        services = self.services

        scale_gs = int(getattr(self, "SCALE_GS", "0"))
        R0_scale = float(getattr(self, "R0_scale", "0"))
        B0_scale = float(getattr(self, "B0_scale", "0"))

        #--- initial force free equilibrium
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
            cur_bc_file = services.get_config_param('CURRENT_BC')
            print('using {}'.format(cur_bc_file))
            instate = Namelist(cur_bc_file)["inbc"]
            if not instate: instate = Namelist(cur_bc_file)["instate"]
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
            shutil.copyfile("g%06d.%05d"%(ishot, itime), cur_eqdsk_file)

        if ps_backend == 'PS':
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':
            ps = plasmastate('ips',1)
            geqdsk = readg(cur_eqdsk_file)
            r0 = geqdsk["rzero" ]
            b0 = abs(geqdsk["bcentr"])
            ip = geqdsk['cpasma']
            ps.init_from_geqdsk (cur_eqdsk_file, nrho=nrho, nth=101)
            inmetric = ps_to_inmetric(ps, r0, b0, ip)

            instate = Namelist(cur_instate_file)
            instate["inmetric"] = inmetric["inmetric"]
            instate.write(cur_instate_file)