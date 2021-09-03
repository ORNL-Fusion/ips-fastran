"""
 -----------------------------------------------------------------------
 elite component
 -----------------------------------------------------------------------
"""

import os
import shutil
from ipsframework import Component
from Namelist import Namelist
from fastran.stability.pdata import pdata

class elite(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('elite.init() entered')

    def step(self, timeid=0):
        print('enter elite.step()')

        #--- code entry
        services = self.services

        ishot = int(services.get_config_param('SHOT_NUMBER'))
        itime = int(services.get_config_param('TIME_ID'))

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- input
        fn_inelite = getattr(self,'INELITE')

        print(fn_inelite)
        inelite = Namelist(fn_inelite)

        nmodes = [int(nmode) for nmode in self.NMODES.split()]

        shutil.copyfile(cur_eqdsk_file, 'eqdsk')

        p = pdata()
        p.load_instate(cur_instate_file)
        p.write('peqdsk')

        #--- run codes
        bin_eq = os.path.join(self.BIN_PATH, self.BIN_EQ)
        bin_vac = os.path.join(self.BIN_PATH, self.BIN_VAC)
        bin_elite = os.path.join(self.BIN_PATH, self.BIN_ELITE)

        for nmode in nmodes:
            inelite['qref_modes']['nn'] = [nmode]

            runid = '%06d_%05d_%02d'%(ishot, itime, nmode)

            inelite.write(runid+".in")

            cwd = services.get_working_dir()

            task_id = services.launch_task(1, cwd, bin_eq+' '+runid, logfile = 'xeq%d.log'%nmode)
            retcode = services.wait_task(task_id)
            if (retcode != 0): raise Exception('Error executing EQ')

            task_id = services.launch_task(1, cwd, bin_vac+' '+runid, logfile = 'xvac%d.log'%nmode)
            retcode = services.wait_task(task_id)
            if (retcode != 0): raise Exception('Error executing VAC')

            task_id = services.launch_task(1, cwd, bin_elite+' '+runid, logfile = 'xelite%d.log'%nmode)
            retcode = services.wait_task(task_id)
            if (retcode != 0): raise Exception('Error executing ELITE')

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('elite.finalize() called')
