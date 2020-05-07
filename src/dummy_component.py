"""
 -----------------------------------------------------------------------
 dummy component
 -----------------------------------------------------------------------
"""

import os
import shutil
from glob import glob
from numpy import *
import netCDF4

#--- ips framework
from component import Component

#--- zcode libraries
from Namelist import Namelist

class dummy_component(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('dummy_component.init() called')

    def step(self, timeid=0):
        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- get plasma state file name

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- workflow
        print('dummy component')

        #--- update plasma state files
        services.update_plasma_state()

        #--- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('dummy_component.finalize() called')
