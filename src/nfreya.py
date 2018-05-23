#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 nfreya component
 -----------------------------------------------------------------------
"""

import os
import subprocess
from numpy import *

from component import Component
from Namelist import Namelist
import nfreya_io

class nfreya(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        services = self.services

        #--- excutable

        nfreya_bin = os.path.join(self.BIN_PATH, self.BIN)
        print nfreya_bin

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file names

        ps_backend = getattr(self,'PS_BACKEND','PS')
        if ps_backend == 'PS':
            cur_state_file = services.get_config_param('CURRENT_STATE')
        elif ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- dakota

        innfreya = Namelist("innfreya",case="lower")

        var_list = [ var.upper() for var in innfreya["innfreya"].keys()]
        for key in var_list:
            for k in range(len( innfreya["innfreya"][key] )):
                try:
                    innfreya["innfreya"][key][k] = float(getattr(self, key+"_%d"%k))
                    print key,k, 'updated'
                except AttributeError:
                    pass
        innfreya.write("innfreya")

        #--- generate input

        f_innfreya = "innfreya"
        try:
            dir_data = self.DIR_DATA
        except:
            dir_data = ''
        print dir_data

        if ps_backend == 'PS':
            nfreya_io.write_inputfiles(cur_state_file, cur_eqdsk_file, f_innfreya,dir_data)
        elif ps_backend == 'INSTATE':
            print 'here'
            nfreya_io.write_inputfiles_instate(cur_instate_file, cur_eqdsk_file, f_innfreya,dir_data)

        #--- run nfreya

        print 'run nfreya'

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, nfreya_bin, logfile = 'xnfreya.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing nfreya')

        #--- get output

        scales = {}
        scales['current'] = float(getattr(self,"SCALE_CURRENT","1.0"))
        scales['particle'] = float(getattr(self,"SCALE_PARTICLE","1.0"))

        if ps_backend == 'PS':
            nfreya_io.update_state(cur_state_file, cur_eqdsk_file,scales)
        elif ps_backend == 'INSTATE':
            nfreya_io.update_instate(cur_instate_file, cur_eqdsk_file,scales)

        #--- update plasma state files

        update_state = int(getattr(self,'UPDATE_STATE','1'))

        if update_state:
            services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def finalize(self, timeStamp=0.0):

        return
