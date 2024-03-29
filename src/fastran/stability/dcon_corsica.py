"""
 -----------------------------------------------------------------------
 dcon component
 -----------------------------------------------------------------------
"""

import sys
import os
import shutil
import re
from ipsframework import Component


class dcon_corsica(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('dcon_corsica.init() called')

    def step(self, timeid=0):
        print('dcon_corsica.step() started')

        # -- excutable
        dcon_bin = os.path.join(self.BIN_PATH, self.BIN)
        print('dcon_bin = ', dcon_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- prepare run
        nmode = int(self.NMODE)
        shutil.copyfile(cur_eqdsk_file, 'eqdsk')

        # -- run dcon
        print('run corsica-dcon')

        try:
            cwd = self.services.get_working_dir()
            task_id = self.services.launch_task(1, cwd, dcon_bin + " s%d %d"%(nmode, nmode), logfile = 'xdcon.log')
            retcode = self.services.wait_task(task_id)
        except Exception:
            print('...in launch_task')
            raise

        if (retcode != 0):
           print('retcode = ', retcode)

        # -- get dcon output
        self.read_dcon("s%d.log"%(nmode))

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('dcon_corsica.finalize() called')

    def read_dcon(self, fn_log, iverb=True):
        for line in open(fn_log,"r").readlines():
            if re.compile("\s*betaN_no-wall").search(line):
                betan_no_wall = float(line.split(":")[-1])
            if re.compile("\s*betaN_ideal-wall").search(line):
                betan_ideal_wall = float(line.split(":")[-1])
        if iverb:
            print("betan_limit = %5.3f %5.3f"%(betan_no_wall, betan_ideal_wall))
        self.betan_no_wall = betan_no_wall
        self.betan_ideal_wall = betan_ideal_wall
