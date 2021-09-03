"""
 -----------------------------------------------------------------------
 cesol driver
 -----------------------------------------------------------------------
"""

import sys
import os
import shutil
from ipsframework import Component

class cesol_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timestamp=0):
        print('>>> cesol_driver.init() called')

    def step(self, timestamp=0):
        print('>>> cesol_driver.step() started')

        #-- entry
        services = self.services

        #-- stage input and plasma state files
        services.stage_input_files(self.INPUT_FILES)
        services.stage_state()

        #-- get list of ports
        ports = services.get_config_param('PORTS')
        port_names = ports['NAMES'].split()
        print('PORTS =', port_names)

        #-- instantiate components in port_names list, except DRIVER itself
        port_dict = {}
        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name)
            port_dict[port_name] = port

        #-- initial time stamp
        t = 0

        #-- initialize components in PORTS list for startup or restart
        init_mode = 'init'
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            services.call(port_dict[port_name], init_mode, t)

        #-- post init processing: stage output
        services.stage_output_files(t, self.OUTPUT_FILES)

        #-- iteration
        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print("number of interation :", nstep)

        for kstep in range(nstep):
            t = kstep
            for port_name in port_names:
                if port_name in ["INIT", "DRIVER"]: continue
                services.call(port_dict[port_name], 'step', t)

        #-- call finalize on each component
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            services.call(port_dict[port_name], 'finalize', t)

    def finalize(self, timestamp = 0):
        print('>>> cesol_driver.finalized() started')

        sym_root = self.services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root, 'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f 0.0\n")
        f.write("%6.3f 0.0\n")
        f.close()
