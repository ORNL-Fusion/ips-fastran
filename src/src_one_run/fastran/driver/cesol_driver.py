"""
 -----------------------------------------------------------------------
 cesol driver
 -----------------------------------------------------------------------
"""

import sys
import os
import shutil
from collections import OrderedDict
from ipsframework import Component

class cesol_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> cesol_driver.init() called')
        print('timeid =', timeid)

        #-- stage input and plasma state files
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.stage_state()

        #-- get list of ports
        ports_names = getattr(self, "PORTS", "").split()
        if len(ports_names) == 0:
            ports_names = self.services.get_config_param('PORTS')['NAMES'].split()
        print('PORTS =', ports_names)

        self.ports = OrderedDict() 
        for port_name in ports_names:
            if port_name in ["INIT", "DRIVER"]: continue
            self.ports[port_name] = self.services.get_port(port_name)

        #-- initialize components in PORTS list 
        init_mode = 'init'
        for port_name in self.ports:
            self.services.call(self.ports[port_name], init_mode, timeid)

        #-- stage output
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def step(self, timeid=0):
        print('>>> cesol_driver.step() started')
        print('timeid =', timeid)

        #-- stage input and plasma state files
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.stage_state()

        #-- iteration
        nstep = int(self.services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print("number of interation :", nstep)

        for kstep in range(nstep):
            print('\n'+72*'='+'\n= cesol driver: iteration number = {}'.format(kstep))

            for port_name in self.ports:
                self.services.call(self.ports[port_name], 'step', kstep)

            self.services.stage_state()
            self.services.stage_output_files(kstep, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('>>> cesol_driver.finalized() started')
        print('timeid =', timeid)

        #-- call finalize on each component
        for port_name in self.ports:
            self.services.call(self.ports[port_name], 'finalize', timeid)

        sym_root = self.services.get_config_param('SIM_ROOT')
        outfile = os.path.join(sym_root, 'RESULT')
        f = open(outfile,"w")
        f.write("0.0\n")
        f.write("0.0\n")
        f.close()
