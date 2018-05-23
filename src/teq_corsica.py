#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 teq component
 -----------------------------------------------------------------------
"""

import os
import shutil
from numpy import *
from component import Component
from plasmastate import plasmastate

class teq(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        pass

    def step(self, timeStamp=0):

        #-- entry

        services = self.services

        #-- excutable

        teq_bin = os.path.join(self.BIN_PATH, self.BIN)
        print 'teq_bin = ',teq_bin

        #-- stage plasma state files

        services.stage_plasma_state()

        #-- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #-- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #-- prepare run

        shutil.copyfile(cur_eqdsk_file,'eqdsk')

        f_inbas = getattr(self,'INBAS')

        #-- run teq

        print 'run corsica-teq'

        try:
            cwd = services.get_working_dir()
            task_id = services.launch_task(1, cwd, teq_bin+" "+f_inbas, logfile = 'xteq.log')
            retcode = services.wait_task(task_id)
        except Exception:
            raise Exception('...in launch_task, teq')

        if (retcode != 0):
            print 'retcode = ',retcode

        #-- get output

        shutil.copy(cur_eqdsk_file+"_inv_teq", cur_eqdsk_file)

        #--- load geqdsk to plasma state file

        ps = plasmastate('ips',1)
        ps.read(cur_state_file)
        ps.load_geqdsk(cur_eqdsk_file)
        ps.store(cur_state_file)

        #-- update plasma state files

        ervices.update_plasma_state()

        #-- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def finalize(self, timeStamp=0):

        return
