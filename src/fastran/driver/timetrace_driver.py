"""
 -----------------------------------------------------------------------
 driver for fastran time-dependent modeling
 -----------------------------------------------------------------------
"""

import os
import glob
import shutil
from ipsframework import Component

class timetrace_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> timetrace_driver.init() started')

        #-- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

    def step(self, timeid=0):
        print('>>> timetrace_driver.driver() started')

        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        #-- get list of ports
        config_port_names = self.services.get_config_param('PORTS')['NAMES']
        port_names = []
        for port_name in config_port_names.split():
            if port_name not in port_names: port_names.append(port_name)

        #-- instantiate components except DRIVER 
        ports = {}
        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = self.services.get_port(port_name)
            if(port == None):
                raise Exception('Error accessing %s component'%port_name)
            ports[port_name] = port

        #-- simulation mode
        sim_mode = self.services.get_config_param('SIMULATION_MODE')
        print('SIMULATION_MODE =', sim_mode)

        #-- init components 
        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            print('init {}'.format(port_name))
            self.services.call(ports[port_name], init_mode, timeid)

        #-- stage plasma state files
        self.services.stage_state()

        #-- stage output after init of each components
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

        #--time loop
        nstep = int(self.services.sim_conf["TIME_LOOP"]["NSTEP"])
        for k in range(timeid, timeid+nstep):
            print('')
            print(72*"=")
            print('= fastran timetrace driver: time step = ', k)
            print('')
            self.services.update_time_stamp(k)

            for port_name in port_names:
                if port_name in ["INIT", "DRIVER"]: continue
                self.services.call(ports[port_name], 'step', k)

            self.services.stage_state()
            self.services.stage_output_files(k, self.OUTPUT_FILES)

        print('')
        print(72*"=")
        print('= fastran driver: finalize')
        print('')
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']: continue
            self.services.call(ports[port_name], 'finalize', k)

    def finalize(self, timeid = 0):
        print(80*'*')
        print('FASTRAN DIRVER FINALIZE')
        sym_root = self.services.get_config_param('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f\n"%0.0)
        f.write("%6.3f\n"%0.0)
        f.close()

        dir_summary = getattr(self, "SUMMARY", "")
        if dir_summary != "":
            if not os.path.exists(dir_summary): os.makedirs(dir_summary)
            dir_state = self.services.get_config_param('STATE_WORK_DIR')
            for filename in glob.glob(os.path.join(dir_state, "*.*")):
                print(filename)
                shutil.copy(filename, dir_summary)
