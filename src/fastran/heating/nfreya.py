"""
 -----------------------------------------------------------------------
 nfreya component
 -----------------------------------------------------------------------
"""

import os
import subprocess
from numpy import *
from Namelist import Namelist
from ipsframework import Component
from fastran.heating import nfreya_io


class nfreya(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('nfreya.init() called')

    def step(self, timeid=0):
        print('nfreya.step() called')

        # -- excutable
        nfreya_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(nfreya_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        if ps_backend == 'PS':
            cur_state_file = self.services.get_config_param('CURRENT_STATE')
            cur_bc_file = self.services.get_config_param('CURRENT_BC')
        elif ps_backend == 'INSTATE':
            cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
            cur_bc_file = cur_instate_file

        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- dakota
        innfreya = Namelist("innfreya", case="lower")
        var_list = [var.upper() for var in innfreya["innfreya"].keys()]
        for key in var_list:
            for k in range(len(innfreya["innfreya"][key])):
                try:
                    innfreya["innfreya"][key][k] = float(getattr(self, key+"_%d" % k))
                    print(key, k, 'updated')
                except AttributeError:
                    pass
        innfreya.write("innfreya")

        FEEDBACK = int(getattr(self, 'FEEDBACK', '10000'))
        inbc = Namelist(cur_bc_file)
        innfreya = Namelist("innfreya", case="lower")
        if timeid >= FEEDBACK:
            pnb = inbc["feedback"]["pnb"][-1]
            innfreya["innfreya"]["bptor"] = [pnb]
            print("feedback:")
            print(inbc["feedback"]["pnb"])
            print(pnb)
        else:
            inbc["feedback"]["pnb"] = [innfreya["innfreya"]["bptor"][0]]
        inbc.write(cur_bc_file)
        innfreya.write("innfreya")

        # -- generate input
        f_innfreya = "innfreya"
        try:
            dir_data = self.DIR_DATA
        except:
            dir_data = ''
        print(dir_data)

        if ps_backend == 'PS':
            nfreya_io.write_inputfiles(cur_state_file, cur_eqdsk_file, f_innfreya, dir_data)
        elif ps_backend == 'INSTATE':
            nfreya_io.write_inputfiles_instate(cur_instate_file, cur_eqdsk_file, f_innfreya, dir_data)

        # -- run nfreya
        print('run nfreya')
        if int(getattr(self, 'SERIAL', '0')) == 1:
            print('nfreya, subprocess')
            logfile = open('xnfreya.log', 'w')
            retcode = subprocess.call([nfreya_bin], stdout=logfile, stderr=logfile, shell=True)
            logfile.close()
        else:
            cwd = self.services.get_working_dir()
            task_id = self.services.launch_task(1, cwd, nfreya_bin, logfile='xnfreya.log')
            retcode = self.services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing nfreya')

        # -- get output
        scales = {}
        scales['current'] = float(getattr(self, "SCALE_CURRENT", "1."))
        scales['particle'] = float(getattr(self, "SCALE_PARTICLE", "1."))

        if ps_backend == 'PS':
            nfreya_io.update_state(cur_state_file, cur_eqdsk_file, scales)
        elif ps_backend == 'INSTATE':
            nfreya_io.update_instate(cur_instate_file, cur_eqdsk_file, scales)

        # -- update plasma state files
        update_state = int(getattr(self, 'UPDATE_STATE', '1'))

        if update_state:
            self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

        # -- clean up
        clean_after = int(getattr(self, "CLEAN_AFTER", "1"))
        if clean_after:
            delete_files = ["bpltfil", "eqpltfil", "namelists", "qikone", "isllog", "runlog", "test_trnspt_mhd.txt"]
            for f in delete_files:
                if os.path.exists(f):
                    os.remove(f)

    def finalize(self, timeid=0):
        print('nfreya.finalize() called')
