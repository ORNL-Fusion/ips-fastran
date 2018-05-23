#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 dummy driver
 -----------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *

from component import Component

class dummy_driver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    # ------------------------------------------------------------------
    # init function

    def init(self, timestamp=0):
        print 80*'*'
        print 'DIRVER INIT'
        return

    # ------------------------------------------------------------------
    # step function

    def step(self, timestamp=0):

        #-- entry

        services = self.services

        #-- stage input and plasma state files

        services.stage_input_files(self.INPUT_FILES)
        services.stage_plasma_state()

        #-- get list of ports

        ports = services.getGlobalConfigParameter('PORTS')
        port_names = ports['NAMES'].split()
        print 'PORTS =', port_names
        port_dict = {}

        print 'HOST =', os.environ["HOST"]

        #-- instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name)
            if(port == None):
                print 'Error accessing %s component'%port_name
                raise
            port_dict[port_name] = port

        #-- initial time stamp

        t = 0

        #-- initialize components in PORTS list for startup or restart

        init_mode = 'init'

        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue
            self.component_call(services,port_name,port_dict[port_name],init_mode,t)

        #-- get plasma state files into driver work directory

        services.stage_plasma_state()

        #-- post init processing: stage output

        services.stage_output_files(t, self.OUTPUT_FILES)

        #-- workflow

        print 'dummy driver called'

        for port_name in port_names:

            if port_name in ["INIT","DRIVER"]: continue
            self.component_call(services,port_name,port_dict[port_name],'step',t)

        #-- call finalize on each component

        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue
            self.component_call(services, port_name, port_dict[port_name], 'finalize', t)

    # ------------------------------------------------------------------
    # finalize function

    def finalize(self, timestamp = 0):

        print 80*'*'
        print 'DIRVER FINALIZE'
        sym_root = self.services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f 0.0\n")
        f.write("%6.3f 0.0\n")
        f.close()

    # ----------------------------------------------------------------------
    # driver methods

    def component_call(self, services, port_name, comp, mode, time):

            comp_mode_string = port_name + ' ' + mode
            print (comp_mode_string)
            try:
                services.call(comp, mode, time)
            except Exception:
                services.exception(comp_mode_string + ' failed')
                raise
            else:
                print comp_mode_string + ' finished'

            return
