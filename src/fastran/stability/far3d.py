"""
 -----------------------------------------------------------------------
 far3d component
 -----------------------------------------------------------------------
"""

import os
import shutil
from Namelist import Namelist
from ipsframework import Component
from fastran.stability import far3d_io
from fastran.util import dakota_io
from fastran.util.fastranutil import freeze


class far3d(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('far3d.init() called')

    def step(self, timeid=0):
        print('far3d.step() started')

        # -- freeze/resume
        if freeze(self, timeid, 'far3d'): return None

        # -- excutable
        far3d_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(far3d_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        f_input_model = 'Input_Model'
        if self.INPUT_MODEL != f_input_model:
            shutil.copyfile(self.INPUT_MODEL, f_input_model)

        f_ineq = self.INPUT_EQ # this is a temporary implementation, will be replaced by a process to generate far3d equilibrium from geqdsk

        # -- dakota binding
        # not implemented

        # -- generate far3d input
        far3d_profile = far3d_io.far3d_io_profile()
        far3d_profile.from_state(f_instate=cur_instate_file, f_state=cur_state_file)
        far3d_profile.write_profile('Profile.txt')

        # -- run genray
        print('run far3d')

        cwd = self.services.get_working_dir()
        task_id = self.services.launch_task(1, cwd, far3d_bin, logfile='xfar3d.log')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing: far3d')

        # placeholder for update state from far3D
        # far3d_io.update_state(cur_state_file, cur_instate_file, cur_eqdsk_file) 

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('genray.finalize() called')
