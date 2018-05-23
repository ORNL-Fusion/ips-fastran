#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 dummy monitor 
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from component import Component

from Namelist import Namelist

class monitor(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):
        return

    def step(self, timeStamp=0.0):
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
