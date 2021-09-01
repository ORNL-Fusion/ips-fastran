"""
 -----------------------------------------------------------------------
 chease component
 -----------------------------------------------------------------------
"""

import os
import shutil
from component import Component
from Namelist import Namelist
from plasmastate import plasmastate

class chease(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('chease.init() started')

        #-- initail equlibrium
        try:
           init_run = int(self.INIT_RUN)
        except:
           init_run = 0
        print('init_run = ',init_run)

        if init_run:
            self.step(-1)

        print('chease.init() done')

    def step(self, timeid=0):
        print('chease.step() started')

        #--- code entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- set default
        NRBOX = int(getattr(self, "NRBOX", "400"))
        NZBOX = int(getattr(self, "NZBOX", "400"))

        #--- input
        chease = Namelist()

        chease["eqdata"]["NEQDSK"] = [1]
        chease["eqdata"]["NSURF" ] = [6]
        chease["eqdata"]["NIDEAL"] = [6]
        chease["eqdata"]["NRBOX" ] = [NRBOX]
        chease["eqdata"]["NZBOX" ] = [NZBOX]

        chease.write("chease_namelist")

        shutil.copyfile(cur_eqdsk_file, "EXPEQ")

        #--- excutables
        chease_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(chease_bin)

        #--- run chease
        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, chease_bin, logfile = 'xchease.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing chease')

        #--- update local geqdsk state
        shutil.copyfile('EQDSK_COCOS_02.OUT', cur_eqdsk_file)

        #--- load geqdsk to plasma state file
        do_update_state = int(getattr(self, "UPDATE_STATE", "0"))

        if do_update_state:
            print('CHEASE: LOAD GEQDSK')
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

        #--- code exit
        print('chease.step() done')

    def finalize(self, timeid=0):
        print('chease.finalize() called')
