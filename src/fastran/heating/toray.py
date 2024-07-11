"""
 -----------------------------------------------------------------------
 toray component
 -----------------------------------------------------------------------
"""

import os
import shutil
import subprocess
import numpy as np
from Namelist import Namelist
from ipsframework import Component
from fastran.heating import toray_io
from fastran.equilibrium.efit_eqdsk import readg
from fastran.plasmastate.plasmastate import plasmastate
from fastran.util import dakota_io
from fastran.state.instate import Instate


class toray(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('toray.init() called')
        return

    def step(self, timeid=0):
        print('toray.step() started')

        # -- stage plasma state files
        self.services.stage_state()

        # -- excutable
        toray_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(toray_bin)

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- get input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- dakota binding
        intoray = Namelist('intoray')

        print('start dakota update')
        dakota_io.update_namelist(self, intoray, section='intoray')
        dakota_io.update_namelist(self, intoray)

        # -- from instate
        if int(getattr(self, 'TRACE', '0')):
            print('update EC power from instate')
            cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
            instate = Instate(cur_instate_file)
            pech = instate['pech']
            intoray['intoray']['rfpow'] = pech

        intoray.write('intoray')

        # -- return if injection power = 0
        # ...

        # -- generate toray input
        geq = readg(cur_eqdsk_file)

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        intoray = Namelist('intoray', case='lower')
        ntoray = intoray['intoray']['ntoray'][0]
        toray_nml = Namelist()
        for key in intoray['edata'].keys():
            toray_nml['edata'][key] = intoray['edata'][key]
        toray_nml.head = 'generated by toray_io.py\n'
        toray_nml.write('toray.in')

        # -- loop over each source of ECH
        cwd = self.services.get_working_dir()

        print('number of gyrotron: ', ntoray)
        for k in range(ntoray):
            toray_io.input_from_state(geq, ps, intoray, k)

            # -- run toray
            if int(getattr(self, 'SERIAL', '0')) == 1:
                print('toray, subprocess')
                logfile = open('xtoray.log', 'w')
                retcode = subprocess.call(
                    [toray_bin], 
                    stdout=logfile, 
                    stderr=logfile, 
                    shell=True)
                logfile.close()
            else:
                task_id = self.services.launch_task(
                    1, 
                    cwd, 
                    toray_bin, 
                    logfile='xtoray.log')
                retcode = self.services.wait_task(task_id)

            if retcode != 0:
                raise Exception('Error executing toray')

            shutil.copyfile('toray.nc', 'toray_%d.nc' % k)

        # -- get toray output
        toray_io.update_state(geq, ps, intoray)

        ps.store(cur_state_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('toray.finalized() called')
