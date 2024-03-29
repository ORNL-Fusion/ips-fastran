"""
 -----------------------------------------------------------------------
 esc component
 -----------------------------------------------------------------------
"""
import os
import shutil
import subprocess
from Namelist import Namelist
from fastran.equilibrium import efit_io
from fastran.plasmastate.plasmastate import plasmastate
from ipsframework import Component

class esc(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('esc.init() entered')
        try:
           init_run = int(self.INIT_RUN)
        except:
           init_run = 0
        print('init_run = ',init_run)

        if init_run:
            self.step(-1)

    def step(self, timeid=0):
        print('enter esc.step()')

        #--- code entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- generate inefit
        efit_io.io_input_from_state(cur_state_file, cur_bc_file)

        #--- excutables
        esc_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(esc_bin)

        wgeqdsk_bin = os.path.join(self.BIN_PATH_WGEQDSK, BIN_WGEQDSK)
        print(wgeqdsk_bin)

        #--- run esc
        cwd = services.get_working_dir()

        if int(getattr(self, 'SERIAL', '0')) == 1:
            print ('esc, subprocess')
            logfile = open('xesc.log', 'w')
            retcode = subprocess.call([esc_bin] ,stdout=logfile, stderr=logfile, shell=True)
            logfile.close()
        else:
            task_id = services.launch_task(1, cwd, esc_bin, logfile = 'xesc.log')
            retcode = services.wait_task(task_id)

        if (retcode != 0):
           print('Error executing ', 'esc')
           raise

        #--- run wgeqdsk
        logfile = open('wgeqdsk.log', 'w')
        retcode = subprocess.call([wgeqdsk_bin], stdout=logfile, stderr=logfile, shell=True)
        logfile.close()

        if (retcode != 0):
           print('Error executing ', 'wgeqdsk')
           raise

        #--- update local geqdsk state
        shutil.copyfile('geqdsk', cur_eqdsk_file)

        #--- load geqdsk to plasma state file
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        ps.load_geqdsk(cur_eqdsk_file)
        ps.store(cur_state_file)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('esc.finalize() called')
