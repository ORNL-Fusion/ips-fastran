"""
 -----------------------------------------------------------------------
 genray component
 -----------------------------------------------------------------------
"""

import os
import shutil
import subprocess
from fastran.heating import genray_io
from Namelist import Namelist
from ipsframework import Component

class genray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('genray.init() called')

        dir_state = self.services.get_config_param('STATE_WORK_DIR')
        open(os.path.join(dir_state, "genray.nc"), "w").close()

    def step(self, timeid=0):
        #--- entry
        print('genray.step() started')
        services = self.services

        ifreeze = int(getattr(self, "FREEZE", -1))
        iresume = int(getattr(self, "RESUME", -1))
        if ifreeze >= 0 and timeid >= ifreeze:
            if iresume < 0 or timeid < iresume:
                print("genray skipped, FREEZE = %d, RESUME = %d, TIMEID = %d"%(ifreeze, iresume, timeid))
                return None

        #--- excutable
        genray_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(genray_bin)

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        if self.INGENRAY != "ingenray":
            shutil.copyfile(self.INGENRAY, "ingenray")

        #--- dakota binding
        f_ingenray = "ingenray"
        ingenray = Namelist(f_ingenray, case="upper")
        for key in ingenray.keys():
            for var in ingenray[key].keys():
                for k in range(len(ingenray[key][var])):
                    if hasattr(self, "%s_%s_%d"%(key, var, k)):
                        ingenray[key][var][k] = float(getattr(self, "%s_%s_%d"%(key, var, k)))
                        print(key, var, k,'updated')
        ingenray.write(f_ingenray)

        #--- generate genray input
        unit = getattr(self, "UNIT", "")
        if unit=="MKS":
            MKS = True
        else:
            MKS = False
        print('GENRAY UNIT:', unit, MKS)

        jmulti = float(getattr(self, "JMULTI", 1.0))
        print('GENRAY JMULTI:', jmulti)

        rho_smooth = float(getattr(self, "RHO_SMOOTH", -1.0))
        print('GENRAY RHO_SMOOTH:', rho_smooth)

        genray_io.write_inputfiles(cur_state_file, cur_eqdsk_file, f_ingenray, MKS)

        add = int(getattr(self, "ADD", "0"))
        print('add = ', add)

        imode = getattr(self, "IMODE", "IC")
        print('imode = ', imode)

        #--- run genray
        print('run genray')

        if int(getattr(self, 'SERIAL','0')) == 1:
            print('genray, subprocess')
            logfile = open('xgenray.log', 'w')
            retcode = subprocess.call([genray_bin], stdout=logfile, stderr=logfile, shell=True)
            logfile.close()
        else:
            cwd = services.get_working_dir()
            task_id = services.launch_task(1, cwd, genray_bin, logfile = 'xgenray.log')
            retcode = services.wait_task(task_id)
        if (retcode != 0):
           raise Exception('Error executing: xgenray')

        #--- get genray output
        genray_io.update_state(cur_state_file, cur_eqdsk_file, imode=imode, jmulti=jmulti, add=add, rho_smooth=rho_smooth)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('genray.finalize() called')
