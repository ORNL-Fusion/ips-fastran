"""
 -----------------------------------------------------------------------
 c1 component
 -----------------------------------------------------------------------
"""

import os
import shutil
import subprocess
import numpy as np
from Namelist import Namelist
from ipsframework import Component
from fastran.util import dakota_io
from fastran.util.fastranutil import freeze


class c1(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('c1.init() called')

    def step(self, timeid=0):
        print('c1.step() started')

        # -- excutable
        c1_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(c1_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- dakota binding
        f_inc1 = self.INC1
        inc1 = Namelist(f_inc1, case='upper')

        print('start dakota update')
        dakota_io.update_namelist(self, inc1)
        inc1.write(f_inc1)

        self.write_inputfile(cur_instate_file, f_inc1)

        # -- run c1
        print('run c1')

        if int(getattr(self, 'SERIAL', '0')) == 1:
            print('c1, subprocess')
            logfile = open('c1.log', 'w')
            retcode = subprocess.call(
                [c1_bin], stdout=logfile, stderr=logfile, shell=True)
            logfile.close()
        else:
            cwd = self.services.get_working_dir()
            task_id = self.services.launch_task(
                1, cwd, c1_bin, logfile='c1.log')
            retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            raise Exception('Error executing: c1')

        # -- get c1 output
        self.update_state(cur_instate_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(
            timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('c1.finalize() called')

    def write_inputfile(self, cur_instate_file, f_inc1, qe_min=0.01, qi_min=0.01):
        instate = Namelist(cur_instate_file)
        inc1 = Namelist(f_inc1)

        qe = instate['instate']['pnbe'][0] \
            + instate['instate']['prfe'][0] \
            + instate['instate']['pfuse'][0] \
            + instate['instate']['pei'][0] \
            + instate['instate']['prad'][0]

        qi = instate['instate']['pnbi'][0] \
            + instate['instate']['prfi'][0] \
            + instate['instate']['pfusi'][0] \
            - instate['instate']['pei'][0]

        area = instate['instate']['area'][0]

        rb = instate["instate"]["rbdry"]
        zb = instate["instate"]["zbdry"]

        i_start = np.argmax(rb)
        i_end = np.argmax(zb)
        print('OMP, X point:', i_start, i_end)

        LX = 0.
        for i in range(i_start, i_end - 1):
            dl = ((rb[i + 1] - rb[i])**2 + (zb[i + 1] - zb[i])**2)**0.5
            LX += dl
        print('LX =', LX)
        print('L0 =', LX / 0.75)

        print('area =', area)
        print('qe =', qe, qe / area)
        print('qi =', qi, qi / area)

        inc1['inc1']['l0'] = [LX / 0.75]
        inc1['inc1']['lx'] = [LX]
        inc1['inc1']['qe'] = [1.e6 * max(qe, qe_min) / area]
        inc1['inc1']['qi'] = [1.e6 * max(qi, qi_min) / area]

        inc1.write(f_inc1)

    def update_state(self, cur_instate_file):
        output = Namelist('sol_summary.dat')

        instate = Namelist(cur_instate_file)
        instate['instate']['np_up'] = output['output']['nu']
        instate['instate']['np_div'] = output['output']['nd']
        instate['instate']['te_div'] = output['output']['te']
        instate['instate']['ti_div'] = output['output']['ti']
        instate['instate']['qe_div'] = [1.e-6 * output['output']['qe_div'][0]]
        instate.write(cur_instate_file)
