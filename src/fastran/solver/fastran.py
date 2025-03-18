"""
 -----------------------------------------------------------------------
 fastran transport solver component
 -----------------------------------------------------------------------
"""

import os
import shutil
import numpy as np
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver import fastran_io_ps
from fastran.solver import fastran_io_instate
from fastran.solver import zdata
from fastran.state.instate import Instate
from fastran.util.fastranutil import freeze


class fastran(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('fastran.init() entered')
        self.icalled = 0

    def step(self, timeid=0):
        print('fastran.step() entered')

        # -- freeze/resume
        if freeze(self, timeid, 'fastran'): return None

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        update_next = getattr(self, 'UPDATE_NEXT', 'disabled')
        print(f'update_next = {update_next}')
        if update_next == 'enabled':
            next_state_file = self.services.get_config_param('NEXT_STATE')
            next_instate_file = self.services.get_config_param('NEXT_INSTATE')

        ps_backend = getattr(self, 'PS_BACKEND', 'pyps').lower()
        update_instate = True if getattr(self, 'UPDATE_INSTATE', 'enabled').lower() == 'enabled' else False
        print('update_instate = ', update_instate)

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- infastran control
        if self.INFASTRAN != 'infastran':
            shutil.copyfile(self.INFASTRAN, 'infastran')

        infastran = Namelist('infastran')
        update_infastran = getattr(self, 'UPDATE_INFASTRAN', '')
        if update_infastran:
            print('UPDATE_INFASTRAN', update_infastran)
            try:
                key, timeid_update, switch_update = update_infastran.split(':')
            except:
                raise Exception('Invalid UPDATE_INFASTRAN format', update_infastran)
            if int(timeid_update) <= int(str(timeid).split('_')[-1]):
                print(f'INFASTARN {key} updated', int(switch_update))
                infastran['infastran'][key] = [int(switch_update)]
        infastran.write('infastran')

        if ps_backend == 'instate':
            fastran_io_instate.write_input(cur_instate_file)
        else:
            fastran_io_ps.write_input(cur_state_file, cur_instate_file, cur_eqdsk_file, f_inprof='inprof')
            if update_next == 'enabled':
                fastran_io_ps.write_input(next_state_file, next_instate_file, cur_eqdsk_file, f_inprof='inprof1')

        # -- run fastran
        fastran_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(fastran_bin)

        ncpu = int(self.NPROC)
        nky = int(self.NPROC_KY)
        n1d = ncpu//nky

        print('ncpu = ', ncpu)
        print('n1d  = ', n1d)
        print('nky  = ', nky)

        cwd = self.services.get_working_dir()
        task_id = self.services.launch_task(ncpu, cwd, fastran_bin, '%d' % n1d, '%d' % nky, logfile='xfastran.log')
        retcode = self.services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing: fastran')

        # -- update local plasma state
        relax = float(getattr(self, 'RELAX', 0.5))
        relax_j = float(getattr(self, 'RELAX_J', 1.))
        fni = float(getattr(self, 'FNI', 1.))
        relax_ip = float(getattr(self, 'RELAX_IP', 0.3))

        if self.icalled >= int(getattr(self, 'ADJUST_IP', 10000)):
            adjust_ip = 1
        else:
            adjust_ip = 0

        if ps_backend == 'instate':
            fastran_io_instate.update_state(
                f_instate=cur_instate_file,
                f_fastran='fastran.nc',
                relax=relax)
        else:
            fastran_io_ps.update_state(
                f_state=cur_state_file,
                f_eqdsk=cur_eqdsk_file,
                f_instate=cur_instate_file,
                f_fastran='fastran.nc',
                update_instate=update_instate,
                relax=relax,
                relax_j=relax_j,
                adjust_ip=adjust_ip,
                fni_target=fni,
                relax_ip=relax_ip)

        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

        self.icalled = self.icalled + 1

    def finalize(self, timeid=0):
        print('fastran.finalize() entered')

        cwd = self.services.get_working_dir()
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(self.services.get_config_param('TIME_ID'))

        dir_state = self.services.get_config_param('STATE_WORK_DIR')
        f = 'f%06d.%05d' % (ishot, itime)
        shutil.copyfile('fastran.nc', os.path.join(dir_state, f))
