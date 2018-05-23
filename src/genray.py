#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 genray component
 -----------------------------------------------------------------------
"""

import os
import shutil
from component import Component
import genray_io
from Namelist import Namelist

class genray(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        return

    def step(self, timeStamp=0):

        #--- entry

        services = self.services

        #--- excutable

        genray_bin = os.path.join(self.BIN_PATH, self.BIN)
        print genray_bin

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        print self.INGENRAY
        try:
            shutil.copyfile(self.INGENRAY,"ingenray")
        except:
            pass

        #--- dakota binding

        f_ingenray = "ingenray"
        ingenray = Namelist(f_ingenray,case="upper")
        for key in ingenray.keys():
            for var in ingenray[key].keys():
                for k in range(len(ingenray[key][var])):
                    if hasattr(self,"%s_%s_%d"%(key,var,k)):
                        ingenray[key][var][k] = float(getattr(self,"%s_%s_%d"%(key,var,k)))
                        print key,var,k,'updated'
        ingenray.write(f_ingenray)

        #--- generate genray input

        try:
            unit = self.UNIT
        except:
            unit = ""
        if unit=="MKS":
            MKS = True
        else:
            MKS = False
        print 'GENRAY UNIT:', unit, MKS

        try:
            jmulti = float(self.JMULTI)
        except:
            jmulti = 1.0
        print 'GENRAY JMULTI:', jmulti

        genray_io.write_inputfiles(cur_state_file, cur_eqdsk_file, f_ingenray, MKS)

        add = int(getattr(self,"ADD","0"))
        print 'add = ',add

        #--- run genray

        print 'run genray'

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, genray_bin, logfile = 'xgenray.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'xgenray'
           raise

        #--- get genray output

        genray_io.update_state(cur_state_file, cur_eqdsk_file, imode='IC', jmulti=jmulti, add=add)

        #--- update plasma state files

        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def finalize(self, timeStamp=0):
        return
