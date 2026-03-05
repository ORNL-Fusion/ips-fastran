"""
 -----------------------------------------------------------------------
 far3d component
 -----------------------------------------------------------------------
"""

import os
import shutil
import glob
import subprocess
from pathlib import Path
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

        f_input_model = 'Input_Model_namelist'
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
        if os.path.exists("woutb"):
            os.remove("woutb")
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
        n_values = range(2, 21)
        template_file = "Input_Model_namelist"
        common_files = ["Profile.txt", "woutb"]
        far3d_bin = os.path.join(self.BIN_PATH,"xfar3d_n")
    
        # Read template once
        with open(template_file, "r") as f:
            template = f.read()
    
        def format_list(name, lst):
            return f"{name} = " + ", ".join(str(x) for x in lst) + ","
        def format_scalar(name, value):
            return f"{name} = {value}"
    
        # Step 1: Create folders and prepare files
        run_dirs = []
        for n in n_values:
            dirname = Path(f"n{n}_case")
            dirname.mkdir(exist_ok=True)
    
            # Copy required files
            for file in common_files:
                shutil.copy(file, dirname / file)
    
            # Generate m/n lists
            q_min = 2.5 if n <= 6 else 3
            m_min = int(round(q_min * n))
            m_max = int(round(4.5 * n))
            m_pos = list(range(m_min, m_max + 1))[:15]
            m_neg = [-m for m in m_pos]
            mmeq = list(range(len(m_pos)))
            mm_vals = m_pos + m_neg + mmeq
            nn_vals = [n]*len(m_pos) + [-n]*len(m_neg) + [0]*len(mmeq)
            nneq = [0]*len(mmeq)
    
            ldim_line = format_scalar("ldim", len(mm_vals))
            leqdim_line = format_scalar("leqdim", len(mmeq))
            customized = template.replace("<MM_LINE>", format_list("mm", mm_vals))\
                                 .replace("<NN_LINE>", format_list("nn", nn_vals))\
                                 .replace("<MMEQ_LINE>", format_list("mmeq", mmeq))\
                                 .replace("<NNEQ_LINE>", format_list("nneq", nneq))\
                                 .replace("<LDIM>", ldim_line)\
                                 .replace("<LEQDIM>", leqdim_line)
    
            with open(dirname / "Input_Model", "w") as f:
                f.write(customized)
    
            run_dirs.append(dirname)
    
        # Step 2: Launch FAR3d in each folder (no wait)
        for d in run_dirs:
            self.services.launch_task(1, str(d), far3d_bin, logfile="xfar3d.log")
    
        print("All runs launched.")
    
        # placeholder for update state from far3D
        # far3d_io.update_state(cur_state_file, cur_instate_file, cur_eqdsk_file) 

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('genray.finalize() called')

