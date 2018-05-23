#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 c2 init component 
 -----------------------------------------------------------------------
"""

import sys,os,os.path,shutil,pickle
import subprocess
from numpy import *
import netCDF4

from component import Component

from Namelist import Namelist

class dummy_init (Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        return

    def step(self, timeStamp):

        services = self.services

        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')
        run_id = services.get_config_param('RUN_ID')

        print 'tokamak_id =', tokamak_id
        print 'shot_number =',shot_number
        print 'run_id =', run_id
    
    def checkpoint(self, timeStamp=0.0):

        print 'dummy_init.checkpoint() called'
        
        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timeStamp, self.RESTART_FILES)
        
    def finalize(self, timeStamp=0.0):

        print 'dummy_init.finalize() called'

