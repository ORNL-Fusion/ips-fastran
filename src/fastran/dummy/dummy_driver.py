"""
 -----------------------------------------------------------------------
 dummy driver
 -----------------------------------------------------------------------
"""
import os
from ipsframework import Component

class dummy_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print(80*'*')
        print('DIRVER INIT')

    def step(self, timeid=0):
        #-- entry
        services = self.services

        #-- stage input and plasma state files
        services.stage_input_files(self.INPUT_FILES)
        services.stage_state()

        #-- get list of ports
        ports = services.get_config_param('PORTS')
        port_names = ports['NAMES'].split()
        print('PORTS =', port_names)
        port_dict = {}

        print('HOST =', os.environ["HOST"])

        #-- instantiate components in port_names list, except DRIVER itself
        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port_dict[port_name] = services.get_port(port_name)

        #-- initial time stamp
        t = 0

        #-- initialize components in PORTS list for startup or restart
        init_mode = 'init'

        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            self.component_call(services, port_name, port_dict[port_name], init_mode, t)

        #-- get plasma state files into driver work directory
        services.stage_state()

        #-- post init processing: stage output
        services.stage_output_files(t, self.OUTPUT_FILES)

        #-- workflow
        print('dummy driver called')

        for port_name in port_names:
            if port_name in ["INIT", "DRIVER"]: continue
            self.component_call(services, port_name, port_dict[port_name], 'step', t)

        #-- call finalize on each component
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            self.component_call(services, port_name, port_dict[port_name], 'finalize', t)

    # ------------------------------------------------------------------
    # finalize function

    def finalize(self, timeid=0):
        print(80*'*')
        print('DIRVER FINALIZE')
        sym_root = self.services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root, 'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f 0.0\n")
        f.write("%6.3f 0.0\n")
        f.close()

    # ----------------------------------------------------------------------
    # driver methods

    def component_call(self, services, port_name, comp, mode, time):
        comp_mode_string = port_name + ' ' + mode
        print(comp_mode_string)
        try:
            services.call(comp, mode, time)
        except Exception:
            services.exception(comp_mode_string + ' failed')
            raise
        else:
            print(comp_mode_string + ' finished')
