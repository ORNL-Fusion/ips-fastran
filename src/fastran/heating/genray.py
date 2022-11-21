"""
 -----------------------------------------------------------------------
 genray component
 -----------------------------------------------------------------------
"""

import os
import shutil
import subprocess
from Namelist import Namelist
from ipsframework import Component
from fastran.heating import genray_io
from fastran.util import dakota_io


class genray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('genray.init() called')

        # -- for CQL3D
        dir_state = self.services.get_config_param('STATE_WORK_DIR')
        open(os.path.join(dir_state, "genray.nc"), "w").close()

    def step(self, timeid=0):
        print('genray.step() started')

        # -- freeze/resume
        ifreeze = int(getattr(self, "FREEZE", -1))
        iresume = int(getattr(self, "RESUME", -1))
        if ifreeze >= 0 and timeid >= ifreeze:
            if iresume < 0 or timeid < iresume:
                print("genray skipped, FREEZE = %d, RESUME = %d, TIMEID = %d" % (ifreeze, iresume, timeid))
                return None

        # -- excutable
        genray_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(genray_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- dakota binding
        f_ingenray = self.INGENRAY
        ingenray = Namelist(f_ingenray, case="upper")

        print("start dakota update")
        dakota_io.update_namelist(self, ingenray)
        ingenray.write(f_ingenray)

        # -- generate genray input
        unit = getattr(self, "UNIT", "")
        if unit == "MKS":
            MKS = True
        else:
            MKS = False
        print('GENRAY UNIT:', unit, MKS)

        jmulti = float(getattr(self, "JMULTI", 1.))
        print('GENRAY JMULTI:', jmulti)

        rho_smooth = float(getattr(self, "RHO_SMOOTH", -1.))
        print('GENRAY RHO_SMOOTH:', rho_smooth)

        genray_io.write_inputfiles(cur_state_file, cur_eqdsk_file, f_ingenray, MKS)

        # -- run genray
        print('run genray')

        if int(getattr(self, 'SERIAL', '0')) == 1:
            print('genray, subprocess')
            logfile = open('xgenray.log', 'w')
            retcode = subprocess.call([genray_bin], stdout=logfile, stderr=logfile, shell=True)
            logfile.close()
        else:
            cwd = self.services.get_working_dir()
            task_id = self.services.launch_task(1, cwd, genray_bin, logfile='xgenray.log')
            retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing: xgenray')

        # -- get genray output
        add = int(getattr(self, "ADD", "0"))
        print('GENRAY ADD  = ', add)

        imode = getattr(self, "IMODE", "IC")
        print('GENRAY IMODE = ', imode)

        genray_io.update_state(cur_state_file, cur_eqdsk_file, imode=imode, jmulti=jmulti, add=add, rho_smooth=rho_smooth)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('genray.finalize() called')
