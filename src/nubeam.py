"""
 -----------------------------------------------------------------------
 nubeam component for steady-state solution
 -----------------------------------------------------------------------
"""

import os
import shutil
from numpy import *
import netCDF4

from  component import Component

from Namelist import Namelist
from plasmastate import plasmastate

class nubeam(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('nubeam.init() called')

        #--- entry
        services = self.services

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- get work directory
        workdir = services.get_working_dir()

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- load nubeam geometry
        print('load innubeam')

        ps = plasmastate('ips',1)
        ps.read(cur_state_file)

        innubeam = Namelist("innubeam")
        var_list = [ var.upper() for var in list(innubeam["nbi_config"].keys())]
        for key in var_list:
            for k in range(len( innubeam["nbi_config"][key] )):
                try:
                    innubeam["nbi_config"][key][k] = float(getattr(self, key+"_%d"%k))
                    print(key,k, 'updated')
                except AttributeError:
                    pass
        innubeam.write("innubeam")

        ps.load_innubeam()

        ps.store(cur_state_file)

        #--- update plasma state files
        services.update_plasma_state()

        #--- generate nubeam input

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
            services.exeception('no ADAS parameter')

        try:
            os.environ['PREACTDIR'] = self.PREACT
        except Exception:
            services.exeception('no PREACT parameter')

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- setup nubeam_comp_exec run
        os.environ['NUBEAM_ACTION'] = 'INIT'
        try:
            del os.environ['FRANTIC_ACTION']
        except:
            pass

        nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(nubeam_bin)

        task_id = services.launch_task(1, workdir, nubeam_bin, logfile = 'log.nubeam_init')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing command:  mpi_nubeam_comp_exec: init ')

    def step(self, timeid=0):
        print('nubeam.step() started')
        ifreeze = int(getattr(self, "FREEZE", -1))
        iresume = int(getattr(self, "RESUME", -1))
        if ifreeze >= 0 and timeid >= ifreeze:
            if iresume < 0 or timeid < iresume:
                print("nubeam skipped, FREEZE = %d, RESUME = %d, TIMEID = %d"%(ifreeze, iresume, timeid))
                return None

        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- get work directory
        workdir = services.get_working_dir()

        #--- set nubeam_comp_exec run
        innubeam = Namelist("innubeam")

        ncpu =  int(self.NPROC)
        nstep = innubeam["nubeam_run"]["nstep"][0]
        navg = innubeam["nubeam_run"]["navg"][0]
        if nstep-navg < 0:
           raise Exception("nubeam.py: nstep < navg")

        dt_nubeam = innubeam["nubeam_run"]["dt_nubeam"][0]

        print("ncpu  = ", ncpu)
        print("dt    = ", dt_nubeam)
        print("nstep = ", nstep)
        print("navg  = ", navg)

        difb_0  = innubeam["nbi_model"]["difb_0"  ][0]
        difb_a  = innubeam["nbi_model"]["difb_a"  ][0]
        difb_in = innubeam["nbi_model"]["difb_in" ][0]
        difb_out= innubeam["nbi_model"]["difb_out"][0]

        try:
            difb_0 = float(getattr(self, "DB"))
            print('DB updated')
        except AttributeError:
            pass
        print('DB: ',difb_0)

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        rho_anom = ps["rho_anom"][:]
        ps["difb_nbi"][:] = difb_a + (difb_0-difb_a)*(1.0-rho_anom**difb_in)**difb_out

        ps.store(cur_state_file)

        shutil.copyfile(cur_state_file, "ps0.nc")

        services.update_plasma_state()

        nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        os.environ['NUBEAM_ACTION'] = 'STEP'
        os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f'%(nstep-navg,dt_nubeam)
        os.environ['STEPFLAG'] = 'TRUE'
        os.environ['NUBEAM_POSTPROC'] = 'summary_test' #'FBM_WRITE'
        try:
            del os.environ['FRANTIC_INIT']
        except:
            pass
        os.environ['FRANTIC_ACTION'] = 'NONE' #'none' #'execute'

        #--- run nstep-navg

        print(os.environ['NUBEAM_REPEAT_COUNT'])
        task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, logfile = 'log.nubeam')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            e = 'Error executing command:  mpi_nubeam_comp_exec: step '
            raise Exception(e)

        if navg > 0:
            print('run navg')

            #--- run navg
            os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f'%(1, dt_nubeam)
            for k in range(navg):
                task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, logfile = 'log.nubeam_%d'%k)
                retcode = services.wait_task(task_id)
                if (retcode != 0):
                    e = 'Error executing command:  mpi_nubeam_comp_exec: step avg '
                    raise Exception(e)
                shutil.copyfile("state_changes.cdf", "state_changes_%d.cdf"%k)

            #--- avgerage
            data = {}
            for k in range(navg):
                filename = "state_changes_%d.cdf"%k
                data[k] = netCDF4.Dataset(filename,'r+',format='NETCDF4')

            vars = []
            for v in list(data[0].variables.keys()):
                if v not in ["ps_partial_update", "version_id"]: vars.append(v)

            for var in vars:
                avg = []
                for k in range(navg):
                    avg.append(data[k].variables[var][:])
                avg = average(array(avg),axis=0)
                data[0].variables[var][:] =  avg

            for k in range(navg):
                data[k].close()

            shutil.copyfile ("state_changes_0.cdf","state_changes.cdf")

        #--- update plasma state
        #shutil.copyfile(cur_state_file, "ps1.nc")

        ps = plasmastate('ips',1)
        ps.read(cur_state_file)
        vol = ps["vol"][:]

        services.merge_current_plasma_state("state_changes.cdf", logfile='log.update_state')

        #---
        services.stage_plasma_state()
        ps = plasmastate('ips',1)
        ps.read(cur_state_file)
        ps["vol"][:] = vol
        ps.update_particle_balance()
        ps.store(cur_state_file)
        services.update_plasma_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

        #shutil.copyfile(cur_state_file,"ps2.nc")

    def finalize(self, timeid=0):
        print('nubeam.finalize() called')
