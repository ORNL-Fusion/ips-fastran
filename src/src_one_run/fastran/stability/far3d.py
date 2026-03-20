"""
 -----------------------------------------------------------------------
 far3d component
 -----------------------------------------------------------------------
"""

import os
import shutil
import glob
import subprocess
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

        # YG--generating woutb equilibrium file for FAR3d
        print(f'Using EQDSK file: {cur_eqdsk_file}')
        # STEP 1: Run vmeclauncher
        print('Running vmeclauncher')
        vmeclauncher_path = os.path.join(self.BIN_PATH, 'vmeclauncher_2.py')
        ret = subprocess.run(['python', vmeclauncher_path, cur_eqdsk_file, '.'],
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in vmeclauncher_2.py')
        
        # STEP 2: Copy VMEC input
        print('Copying input.vmec file')
        vmec_input_files = glob.glob('./output/input*')
        if not vmec_input_files:
            raise Exception('No input.vmec file found in ./output/')
        shutil.copy(vmec_input_files[0], 'input.vmec')
        
        # STEP 3: Run xvmec
        print('Running VMEC...')
        ret = subprocess.run(['srun', '-n', '16', 'xvmec', 'input.vmec'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in xvmec')
        
        # STEP 4: Run Booz_xform
        print('Running Booz_xform...')
        ret = subprocess.run(['srun', 'xbooz_xform', 'output_file.txt', 'far'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in xbooz_xform')
        
        print('All preprocessing steps completed successfully.')      
        #f_ineq = self.INPUT_EQ # this is a temporary implementation, will be replaced by a process to generate far3d equilibrium from geqdsk


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

