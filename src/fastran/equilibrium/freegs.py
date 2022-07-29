"""
 -----------------------------------------------------------------------
 freegs component
 -----------------------------------------------------------------------
"""
import os
import shutil
import time as timer
import numpy as np
from Namelist import Namelist
from fastran.equilibrium import efit_io
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver.inmetric_io import ps_to_inmetric
from fastran.equilibrium.efit_eqdsk import readg
from ipsframework import Component

from fastran.equilibrium.freegs_io import call_freegs

class freegs(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def _get_shot_time(self):
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(float(self.TIME_ID)) if hasattr(self, "TIME_ID") else int(self.services.get_config_param('TIME_ID'))
        return ishot, itime 

    def init(self, timeid=0):
        print('>>> freegs.init() started')

        #--- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid
        print('freegs time = {}'.format(itime))

        #--- initial equilibrium
        init_run = int(getattr(self, "INIT_RUN", 0))
        print('init_run = ', init_run)
        if init_run:
            print('initial approximate equilibrium')
            self.services.stage_state()
            self.initial_equilibrium(ishot, itime)
            self.services.update_state()

    def step(self, timeid=0):
        print('>>> freegs.step() started')

        #--- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid
        print('freegs time = {}'.format(itime))

        #--- freeze
        ifreeze = int(getattr(self, "FREEZE", 10000))
        iresume = int(getattr(self, "RESUME", -10000))
        if timeid > ifreeze and timeid < iresume:
            print("FREEZE: timeid = {}, ifreeze = {}, iresume = {}"%(timeid, ifreeze, iresume))
            return

        #--- stage plasma state files
        self.services.stage_state()

        #--- get plasma state file names
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        if ps_backend == 'PS':
            cur_state_file = self.services.get_config_param('CURRENT_STATE')
            cur_bc_file = self.services.get_config_param('CURRENT_BC')
        #elif ps_backend == 'INSTATE':
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        use_instate = int(getattr(self, "USE_INSTATE", "0"))

        #--- generate inefit
        f_inefit = "inefit"
        mode = getattr(self, "PRESSURE", "kinetic")
        print('mode =', mode)

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
                inefit['inefit'][key] = 0.5*np.array(inefit0['inefit'][key]) + 0.5*np.array(inefit['inefit'][key])
            inefit.write(f_inefit)

        #--- topology
        topology = getattr(self, "TOPOLOGY", "")
        print('topology =', topology)

        #--- initial equilibrium for inner iteration
        init_run_step = int(getattr(self, "INIT_RUN_STEP", 0))
        print('init_run_step =', init_run_step)

        cwd = self.services.get_working_dir()

        error = float(getattr(self, "ERROR", "1.0e-4"))
        print('error =', error)

        niter = int(getattr(self, "NITER", "5"))
        print('niter = ', niter)

        iconv = 0
        for k in range(niter):
            print ("*** generate kfile %d"%k)
            inefit = Namelist(f_inefit, "r")

            t0 = timer.time()

            relax = 1 if k > 0 else 0

            ps = plasmastate('ips',1)

            geqdsk = readg(cur_eqdsk_file)

            #r0  = geqdsk["rzero" ]
            #b0  = abs(geqdsk["bcentr"])
            #ip  = geqdsk['cpasma']
            r0 = inefit["inefit"]["r0"][0]
            b0 = inefit["inefit"]["b0"][0]
            ip = inefit["inefit"]["ip"][0]
            print ('r0 = ',r0)
            print ('b0 = ',b0)
            print ('ip = ',ip)

            try:
                ps.init_from_geqdsk(cur_eqdsk_file, nrho=nrho, nth=101)
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
                    inmetric["inmetric"][key] = 0.5*np.array(inmetric_prev["inmetric"][key]) + 0.5*np.array(inmetric["inmetric"][key])

            iconv = efit_io.fixbdry_kfile(ishot, itime, Namelist(f_inefit, "r"), inmetric["inmetric"], relax=relax, topology=topology, error=error)

            if iconv:
                print ('converged')
                break

            t1 = timer.time()
            print("run freegs")
            call_freegs(cur_instate_file, f_inefit, init=False)
            t2 = timer.time()
            print('elapsed time = %6.3f %6.3f'%(t1-t0, t2-t1))

            shutil.copyfile("lsn.geqdsk", cur_eqdsk_file)

            #--- update local geqdsk state
            shutil.copyfile("lsn.geqdsk", "g%06d.%05d"%(ishot, itime))

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
        self.services.update_state()

        #--- archive output files
        if self.OUTPUT_FILES:
            self.services.stage_output_files(timeid, self.OUTPUT_FILES)
        else:
            output_files = ' '.join(["{0}{1:06d}.{2:05d}".format(header, ishot, itime) for header in ["g", "a", "k", "m"]])
            output_files += " inefit"
            print("stage_output_files: ", output_files)
            self.services.stage_output_files(timeid,  output_files)

    def finalize(self, timeid):
        print ('>>> freegs.finalize() called')

    def initial_equilibrium(self, ishot, itime):
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        efit_io.io_input_from_instate(f_instate=cur_instate_file, f_inefit="inefit", mode='pmhd')

        call_freegs(cur_instate_file, "inefit", init=True)

        shutil.copyfile("lsn.geqdsk", cur_eqdsk_file)

        ps = plasmastate('ips',1)
        ps.read(cur_state_file)
        ps.load_geqdsk(cur_eqdsk_file)
        ps.store(cur_state_file)

