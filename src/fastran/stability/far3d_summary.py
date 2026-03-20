"""
 -----------------------------------------------------------------------
 far3d component
 -----------------------------------------------------------------------
"""

import os
import glob
import shutil
import subprocess
from pathlib import Path

from Namelist import Namelist  # if unused, you can remove this import
from ipsframework import Component
from fastran.stability import far3d_io
from fastran.util import dakota_io  # if unused, you can remove this import
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
        if freeze(self, timeid, 'far3d'):
            return None

        # Resolve IPS working directory (absolute)
        cwd = Path(self.services.get_working_dir()).resolve()
        print(f"IPS working dir: {cwd}")

        # -- executable paths
        bin_path = Path(self.BIN_PATH).resolve()
        default_exe = (bin_path / self.BIN).resolve()
        print(f"default exe (from config BIN): {default_exe}")

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # Ensure template present locally
        f_input_model = 'Input_Model_namelist'
        if self.INPUT_MODEL != f_input_model:
            shutil.copyfile(self.INPUT_MODEL, f_input_model)

        # ===== Preprocessing to produce woutb (or equivalent) =====
        print(f'Using EQDSK file: {cur_eqdsk_file}')

        # STEP 1: Run vmeclauncher
#        print('Running vmeclauncher')
        vmeclauncher_path = (bin_path / 'vmeclauncher_2.py').resolve()
        python3 vmeclauncher.path cur_eqdsk_file .
#        ret = subprocess.run(
#            ['python', str(vmeclauncher_path), cur_eqdsk_file, '.'],
#            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
#        )
#        print(ret.stdout)
#        if ret.returncode != 0:
#            print(ret.stderr)
#            raise Exception('Error in vmeclauncher_2.py')

        # STEP 2: Copy VMEC input
        print('Copying input.vmec file')
        vmec_input_files = glob.glob('./output/input*')
        if not vmec_input_files:
            raise Exception('No input.vmec file found in ./output/')
        shutil.copy(vmec_input_files[0], 'input.vmec')

        # STEP 3: Run xvmec
        print('Running VMEC...')
        ret = subprocess.run(
            ['srun', '-n', '16', 'xvmec', 'input.vmec'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in xvmec')

        # ensure any old woutb is removed so downstream steps are fresh
        if Path("woutb").exists():
            Path("woutb").unlink()

        # STEP 4: Run Booz_xform
        print('Running Booz_xform...')
        ret = subprocess.run(
            ['srun', 'xbooz_xform', 'output_file.txt', 'far'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in xbooz_xform')

        print('All preprocessing steps completed successfully.')

        # ===== Generate FAR3d profile input =====
        far3d_profile = far3d_io.far3d_io_profile()
        far3d_profile.from_state(f_instate=cur_instate_file, f_state=cur_state_file)
        far3d_profile.write_profile('Profile.txt')

        # ===== Prepare and launch FAR3d scans =====
        print('run far3d')
        n_values = range(2, 21)
        template_file = "Input_Model_namelist"
        common_files = ["Profile.txt", "woutb"]
        far3d_exe = (bin_path / "xfar3d_n").resolve()

        # Read template once
        with open(template_file, "r") as f:
            template = f.read()

        def format_list(name, lst):
            return f"{name} = " + ", ".join(str(x) for x in lst) + ","

        def format_scalar(name, value):
            return f"{name} = {value}"

        # Step 1: Create folders and prepare files (ABSOLUTE dirs)
        run_dirs = []
        for n in n_values:
            dirname = (cwd / f"n{n}_case")
            dirname.mkdir(parents=True, exist_ok=True)

            # Copy required files (from current CWD into run dir)
            for file in common_files:
                shutil.copy(file, dirname / file)

            # Generate m/n lists
            q_min = 2.5 # if n <= 6 else 3
            m_min = int(round(q_min * n))
            m_max = int(round(3.8 * n))
            m_pos = list(range(m_min, m_max + 1))[:13]
            m_neg = [-m for m in m_pos]
            mmeq = list(range(len(m_pos)))
            mm_vals = m_pos + m_neg + mmeq
            nn_vals = [n] * len(m_pos) + [-n] * len(m_neg) + [0] * len(mmeq)
            nneq = [0] * len(mmeq)

            ldim_line = format_scalar("ldim", len(mm_vals))
            leqdim_line = format_scalar("leqdim", len(mmeq))
            customized = (
                template.replace("<MM_LINE>", format_list("mm", mm_vals))
                        .replace("<NN_LINE>", format_list("nn", nn_vals))
                        .replace("<MMEQ_LINE>", format_list("mmeq", mmeq))
                        .replace("<NNEQ_LINE>", format_list("nneq", nneq))
                        .replace("<LDIM>", ldim_line)
                        .replace("<LEQDIM>", leqdim_line)
            )

            with open(dirname / "Input_Model", "w") as f:
                f.write(customized)

            run_dirs.append(dirname)

        # Number of Dask workers/nodes (from config if available)
        try:
            dask_nodes = int(self.services.get_config_param("DASK_NODES"))
        except Exception:
            dask_nodes = 1
        print(f"DASK_NODES = {dask_nodes}")

        # Step 2: Launch FAR3d across all run directories and wait for completion
        self.run_far3d_with_dask(run_dirs, str(far3d_exe), dask_nodes)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def run_far3d_with_dask(self, run_dirs, far3d_exe, dask_nodes):
        """Enqueue one FAR3d task per run directory, submit via Dask, and wait."""
        if dask_nodes is None or dask_nodes < 1:
            raise Exception("DASK_NODES undefined or < 1")

        pool_name = "far3d_pool"
        self.services.create_task_pool(pool_name)

        # Belt-and-suspenders: assert all dirs exist (absolute)
        for d in run_dirs:
            if not Path(d).is_dir():
                raise Exception(f"Run dir missing before submit: {d}")

        # Enqueue one task per run directory (ABSOLUTE working dir + logfile)
        for i, rdir in enumerate(run_dirs, start=1):
            rdir_abs = str(Path(rdir).resolve())
            logfile_abs = str((Path(rdir) / f"xfar3d_{i:02d}.log").resolve())

            self.services.add_task(
                pool_name,
                f"far3d_task_{i:02d}",  # task label
                1,                      # processes per task
                rdir_abs,               # ABSOLUTE working directory
                far3d_exe,              # executable (absolute path or on PATH)
                logfile=logfile_abs     # ABSOLUTE logfile path
                # , env=my_env_dict     # optional: per-task environment
            )

        # Submit the whole pool to Dask workers
        ret_val = self.services.submit_tasks(
            pool_name,
            use_dask=True,
            dask_nodes=dask_nodes
        )
        print("submit_tasks ret_val =", ret_val)

        # Block until all tasks complete; returns {task_label: exit_code}
        exit_status = self.services.get_finished_tasks(pool_name)
        print("exit_status =", exit_status)

        # Check for failures
        failed = {k: v for k, v in exit_status.items() if v != 0}
        if failed:
            details = ", ".join(f"{k}:{v}" for k, v in failed.items())
            raise Exception(f"FAR3d failures in Dask pool: {details}")

        print("All runs completed.")

    def finalize(self, timeid=0):
        print('far3d.finalize() called')
                
        # Block until all tasks complete; returns {task_label: exit_code}
        exit_status = self.services.get_finished_tasks(pool_name)
        print("exit_status =", exit_status)
        
        # Check for failures
        failed = {k: v for k, v in exit_status.items() if v != 0}
        if failed:
            details = ", ".join(f"{k}:{v}" for k, v in failed.items())
            raise Exception(f"FAR3d failures in Dask pool: {details}")
            
        print("All runs completed.")
        
    def finalize(self, timeid=0):
        print('far3d.finalize() called')
        
        collect_path = (Path(self.BIN_PATH) / 'collect.py').resolve()
        ret = subprocess.run(
            ['python3', str(collect_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        print(ret.stdout)
        if ret.returncode != 0:
            print(ret.stderr)
            raise Exception('Error in collect.py')
 
