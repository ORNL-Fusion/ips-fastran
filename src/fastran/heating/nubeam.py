"""
 -----------------------------------------------------------------------
 nubeam component for steady-state solution
 -----------------------------------------------------------------------
"""

import os
import shutil
from numpy import *
import netCDF4
from Namelist import Namelist
from fastran.plasmastate.plasmastate import plasmastate
from ipsframework import Component
from fastran.util import dakota_io
from fastran.state.instate import Instate


class nubeam(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('nubeam.init() called')

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- get work directory
        workdir = self.services.get_working_dir()

        # -- stage plasma state files
        self.services.stage_state()

        load_init_eq = getattr(self, "LOAD_INIT_EQ", "disabled")
        if load_init_eq == "enabled":
            print("load equilibrium")
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        # -- dakota binding
        innubeam = Namelist("innubeam")

        print("start dakota update")
        dakota_io.update_namelist(self, innubeam, section="nbi_config")
        dakota_io.update_namelist(self, innubeam)

        innubeam.write("innubeam")

        # -- load nubeam geometry
        print('load innubeam')

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        ps.load_innubeam()

        ps.update_particle_balance()
        ps.store(cur_state_file)

        # -- update plasma state files
        self.services.update_state()

        # -- generate nubeam input

        nubeam_files = Namelist()
        nubeam_files["NUBEAM_FILES"]["INPUT_PLASMA_STATE"] = [cur_state_file]
        nubeam_files["NUBEAM_FILES"]["PLASMA_STATE_UPDATE"] = ["state_changes.cdf"]
        nubeam_files["NUBEAM_FILES"]["INIT_NAMELIST"] = ["nubeam_init_input.dat"]
        nubeam_files.write("nubeam_init_files.dat")

        nubeam_files = Namelist()
        nubeam_files["NUBEAM_FILES"]["INPUT_PLASMA_STATE"] = [cur_state_file]
        nubeam_files["NUBEAM_FILES"]["PLASMA_STATE_UPDATE"] = ["state_changes.cdf"]
        nubeam_files["NUBEAM_FILES"]["STEP_NAMELIST"] = ["nubeam_step_input.dat"]
        nubeam_files.write("nubeam_step_files.dat")

        innubeam = Namelist("innubeam")

        nubeam_init_input = Namelist()
        nubeam_init_input["NBI_INIT"] = innubeam["NBI_INIT"]
        nubeam_init_input.write("nubeam_init_input.dat")

        nubeam_step_input = Namelist()
        nubeam_step_input["NBI_UPDATE"] = innubeam["NBI_UPDATE"]
        nubeam_step_input.write("nubeam_step_input.dat")

        try:
            os.environ['ADASDIR'] = self.ADAS
        except Exception:
            self.services.exeception('no ADAS parameter')

        try:
            os.environ['PREACTDIR'] = self.PREACT
        except Exception:
            self.services.exeception('no PREACT parameter')

        # -- stage plasma state files
        self.services.stage_state()

        # -- setup nubeam_comp_exec run
        os.environ['NUBEAM_ACTION'] = 'INIT'
        try:
            del os.environ['FRANTIC_ACTION']
        except:
            pass

        nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(nubeam_bin)

        task_id = self.services.launch_task(1, workdir, nubeam_bin, logfile='log.nubeam_init')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing command:  mpi_nubeam_comp_exec: init ')

    def step(self, timeid=0):
        print('nubeam.step() started')
        ifreeze = int(getattr(self, "FREEZE", -1))
        iresume = int(getattr(self, "RESUME", -1))
        if ifreeze >= 0 and timeid >= ifreeze:
            if iresume < 0 or timeid < iresume:
                print("nubeam skipped, FREEZE = %d, RESUME = %d, TIMEID = %d" % (ifreeze, iresume, timeid))
                return None

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        # self.services.stage_input_files(self.INPUT_FILES)

        # -- get work directory
        workdir = self.services.get_working_dir()

        # -- set nubeam_comp_exec run
        innubeam = Namelist("innubeam")

        ncpu = int(self.NPROC)
        nstep = innubeam["nubeam_run"]["nstep"][0]
        navg = innubeam["nubeam_run"]["navg"][0]
        if nstep-navg < 0:
            raise Exception("nubeam.py: nstep < navg")

        dt_nubeam = innubeam["nubeam_run"]["dt_nubeam"][0]

        print("ncpu  = ", ncpu)
        print("dt    = ", dt_nubeam)
        print("nstep = ", nstep)
        print("navg  = ", navg)

        difb_0 = innubeam["nbi_model"]["difb_0"][0]
        difb_a = innubeam["nbi_model"]["difb_a"][0]
        difb_in = innubeam["nbi_model"]["difb_in"][0]
        difb_out = innubeam["nbi_model"]["difb_out"][0]

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        rho_anom = ps["rho_anom"][:]
        ps["difb_nbi"][:] = difb_a + (difb_0 - difb_a)*(1. - rho_anom**difb_in)**difb_out

        ps.update_particle_balance()  # <--------

        if getattr(self, "TIMEBC", "") == "INSTATE":
           print("innubeam updated from instate")
           instate = Instate(cur_instate_file)
           nbeam = innubeam['NBI_CONFIG']['NBEAM'][0]
           for k in range(nbeam):
              pnbi_k = instate['PNBI_%d'%k][0]
              print(k, pnbi_k)
              innubeam['NBI_CONFIG']['PINJA'][k] = pnbi_k
              ps.load_innubeam()

        ps.store(cur_state_file)

        self.services.update_state()

        nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        os.environ['NUBEAM_ACTION'] = 'STEP'
        os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f' % (nstep-navg, dt_nubeam)
        os.environ['STEPFLAG'] = 'TRUE'
        os.environ['NUBEAM_POSTPROC'] = 'summary_test'  # 'FBM_WRITE'
        try:
            del os.environ['FRANTIC_INIT']
        except:
            pass
        os.environ['FRANTIC_ACTION'] = getattr(self, "FRANTIC_ACTION", "NONE")

        # -- run nstep - navg
        print(os.environ['NUBEAM_REPEAT_COUNT'])
        task_id = self.services.launch_task(self.NPROC, workdir, nubeam_bin, logfile='log.nubeam')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            e = 'Error executing command:  mpi_nubeam_comp_exec: step '
            raise Exception(e)

        # -- run navg
        if navg > 0:
            print('run navg')

            # --- run navg
            os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f' % (1, dt_nubeam)
            for k in range(navg):
                task_id = self.services.launch_task(self.NPROC, workdir, nubeam_bin, logfile='log.nubeam_%d' % k)
                retcode = self.services.wait_task(task_id)
                if (retcode != 0):
                    e = 'Error executing command:  mpi_nubeam_comp_exec: step avg '
                    raise Exception(e)
                shutil.copyfile("state_changes.cdf", "state_changes_%d.cdf" % k)

            #--- avgerage
            data = {}
            for k in range(navg):
                filename = "state_changes_%d.cdf" % k
                data[k] = netCDF4.Dataset(filename, 'r+', format='NETCDF4')

            vars = []
            for v in list(data[0].variables.keys()):
                if v not in ["ps_partial_update", "version_id"]:
                    vars.append(v)

            for var in vars:
                avg = []
                for k in range(navg):
                    avg.append(data[k].variables[var][:])
                avg = average(array(avg), axis=0)
                data[0].variables[var][:] = avg

            for k in range(navg):
                data[k].close()

            shutil.copyfile("state_changes_0.cdf", "state_changes.cdf")

        # -- update plasma state
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        vol = ps["vol"][:]
        self.services.merge_current_state("state_changes.cdf", logfile='log.update_state')

        # --
        self.services.stage_state()
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        ps["vol"][:] = vol
        ps.update_particle_balance()
        ps.store(cur_state_file)
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('nubeam.finalize() called')
