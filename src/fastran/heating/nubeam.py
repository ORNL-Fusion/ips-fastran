"""
 -----------------------------------------------------------------------
 nubeam component for steady-state solution
 -----------------------------------------------------------------------
"""
import os
import shutil
import numpy as np
import netCDF4
from Namelist import Namelist
from ipsframework import Component
from fastran.heating import nubeam_io
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate
from fastran.util import dakota_io
from fastran.util.fastranutil import freeze
from fastran.util.input_default import input_default


class nubeam(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>> nubeam.init() started')

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- get work directory
        workdir = self.services.get_working_dir()

        # -- stage plasma state files
        self.services.stage_state()

        load_init_eq = getattr(self, 'LOAD_INIT_EQ', 'disabled')
        if load_init_eq == 'enabled':
            print('load equilibrium')
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        # -- dakota binding
        innubeam = Namelist('innubeam')

        print('start dakota update')
        dakota_io.update_namelist(self, innubeam, section='nbi_config')
        dakota_io.update_namelist(self, innubeam)

        # -- from instate
        #if int(getattr(self, 'TRACE', '0')):
        #    print('update NB:workdir power from instate')
        #    cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        #    instate = Instate(cur_instate_file)
        #    instate.to_ps_pnbi(ps) 

        innubeam.write('innubeam')

        # -- load to plasma state
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        print('load innubeam')
        ps.load_innubeam()

        print('load difb')
        difb_0 = innubeam['nbi_model']['difb_0'][0]
        difb_a = innubeam['nbi_model']['difb_a'][0]
        difb_in = innubeam['nbi_model']['difb_in'][0]
        difb_out = innubeam['nbi_model']['difb_out'][0]

        rho_anom = ps['rho_anom'][:]
        ps['difb_nbi'][:] = difb_a + (difb_0 - difb_a) * (1. - rho_anom**difb_in)**difb_out

        ps.update_particle_balance()
        ps.store(cur_state_file)

        # -- update plasma state files
        self.services.update_state()

        # -- generate nubeam input

        nubeam_files = Namelist()
        nubeam_files['NUBEAM_FILES']['INPUT_PLASMA_STATE'] = [cur_state_file]
        nubeam_files['NUBEAM_FILES']['PLASMA_STATE_UPDATE'] = ['state_changes.cdf']
        nubeam_files['NUBEAM_FILES']['INIT_NAMELIST'] = ['nubeam_init_input.dat']
        nubeam_files.write('nubeam_init_files.dat')

        nubeam_files = Namelist()
        nubeam_files['NUBEAM_FILES']['INPUT_PLASMA_STATE'] = [cur_state_file]
        nubeam_files['NUBEAM_FILES']['PLASMA_STATE_UPDATE'] = ['state_changes.cdf']
        nubeam_files['NUBEAM_FILES']['STEP_NAMELIST'] = ['nubeam_step_input.dat']
        nubeam_files.write('nubeam_step_files.dat')

        nubeam_init_input = Namelist()
        nubeam_init_input['NBI_INIT'] = innubeam['NBI_INIT']
        nubeam_init_input.write('nubeam_init_input.dat')

        nubeam_step_input = Namelist()
        nubeam_step_input['NBI_UPDATE'] = innubeam['NBI_UPDATE']
        nubeam_step_input.write('nubeam_step_input.dat')

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
        print('>> nubeam.step() started')

        # -- freeze/resume
        if freeze(self, timeid, 'nubeam'): return None

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- get work directory
        workdir = self.services.get_working_dir()

        # -- set nubeam_comp_exec run
        innubeam = Namelist('innubeam')

        ncpu = int(self.NPROC)
        nstep = innubeam['nubeam_run']['nstep'][0]
        navg = innubeam['nubeam_run']['navg'][0]
        if nstep-navg < 0:
            raise Exception('nubeam.py: nstep < navg')

        dt_nubeam = innubeam['nubeam_run']['dt_nubeam'][0]

        print('ncpu  = ', ncpu)
        print('dt    = ', dt_nubeam)
        print('nstep = ', nstep)
        print('navg  = ', navg)

        # -- enforce particle balance
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        ps.update_particle_balance()
        ps.store(cur_state_file)
        self.services.update_state() # needed since the merge_current_state will be applied to STATE_WORK_DIR/CURRENT_STATE

        nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        os.environ['NUBEAM_ACTION'] = 'STEP'
        os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f' % (nstep - navg, dt_nubeam)
        os.environ['STEPFLAG'] = 'TRUE'
        os.environ['NUBEAM_POSTPROC'] = 'summary_test'  # 'FBM_WRITE'
        try:
            del os.environ['FRANTIC_INIT']
        except:
            pass
        os.environ['FRANTIC_ACTION'] = getattr(self, 'FRANTIC_ACTION', 'NONE')

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
                    e = 'Error executing command:  mpi_nubeam_comp_exec: step avg'
                    raise Exception(e)
                shutil.copyfile('state_changes.cdf', 'state_changes_%d.cdf' % k)

            #--- avgerage
            data = {}
            for k in range(navg):
                filename = 'state_changes_%d.cdf' % k
                data[k] = netCDF4.Dataset(filename, 'r+', format='NETCDF4')

            vars = []
            for v in list(data[0].variables.keys()):
                if v not in ['ps_partial_update', 'version_id']:
                    vars.append(v)

            for var in vars:
                avg = []
                for k in range(navg):
                    avg.append(data[k].variables[var][:])
                avg = np.average(np.array(avg), axis=0)
                data[0].variables[var][:] = avg

            for k in range(navg):
                data[k].close()

            shutil.copyfile('state_changes_0.cdf', 'state_changes.cdf')


        # -- update plasma state
        self.services.merge_current_state('state_changes.cdf', logfile='log.update_state')

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        ps.update_particle_balance()
        ps.store(cur_state_file)

        # -- dump to instate
        update_instate = input_default(self, key='UPDATE_INSTATE', default='disabled', alias=[], maps={'0':'disabled', '1':'enabled'})
        if update_instate:
           nubeam_io.update_instate(cur_state_file, cur_instate_file) 

        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('nubeam.finalize() called')
