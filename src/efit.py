#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran efit component 
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

from  component import Component

#-------------------
#--- zcode libraries
import zefit

class efit(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0.0):

        return

    def step(self, timeStamp=0.0):

        #--- entry

        print 'enter efit.step()'

        if (self.services == None) :
            print 'Error in efit.step () : No services'
            raise Exception('Error in efit.step (): No services')
        services = self.services

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
        except Exception, e:
            print 'Error in call to stage_plasma_state()', e

        #--- get plasma state file names

        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk = services.get_config_param('CURRENT_EQDSK')

        #--- generate inefit

        zefit.io_input_from_instate(cur_instate_file)

        #--- run efit

        try:
            efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            efit_bin = os.path.join(self.BIN, 'efitd90 129 129')
        print fastran_bin

        niter = 5
        shot=0
        time=0
        f_inefit="inefit"

        logfile = open('efit.log', 'w')
        logfile.close()
        for k in range(niter):

            print "generate kfile"
            if k == 0:
                zefit.fixbdry_kfile_init(shot,time,f_inefit)
            else:
                zefit.fixbdry_kfile(shot,time,f_inefit)
    
            print "run efit"
            kfile = "k000000.00000"
            args = "2\n 1\n "+kfile
            command = 'echo \"%s\"'%args + ' | ' + efit_bin 

            logfile = open('efit.log', 'a')
            retcode = subprocess.call([command]
                          ,stdout=logfile,stderr=logfile,shell=True)
            if (retcode != 0):
               print 'Error executing ', 'efit'
               raise
            logfile.close()

        #--- update local geqdsk state

        shutil.copyfile('g000000.00000', cur_eqdsk)

        #--- update plasma state files

        try:
            services.update_plasma_state()
        except Exception, e:
            print 'Error in call to update_plasma_state()', e
            raise

        #--- archive output files

        try:
            services.stage_output_files(timeStamp, self.OUTPUT_FILES)
        except Exception, e:
            print 'Error in call to stage_output_files()', e
            raise Exception, e
 
        return
    
    def finalize(self, timeStamp=0.0):
        return
    
