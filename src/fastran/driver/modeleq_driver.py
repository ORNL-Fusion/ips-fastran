"""
 fastran model equilibrium driver
"""

import os
import glob
import shutil
import time as timer
from Namelist import Namelist
from ipsframework import Component


class modeleq_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('* MODELEQ DIRVER INIT')

    def step(self, timeid=0):
        # -- stage input and plasma state files
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.stage_state()

        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        # -- get list of ports
        ports = self.services.get_config_param('PORTS')
        port_names = ports['NAMES'].split()
        print('PORTS =', port_names)
        port_dict = {}

        post_names = ports['POSTS'].split()
        print(('POSTS =', post_names))

        # -- instantiate components in port_names list, except DRIVER itself
        for port_name in port_names:
            if port_name in ["DRIVER"]:
                continue
            port_dict[port_name] = self.services.get_port(port_name)

        # -- initialize components in PORTS list for startup
        init_mode = 'init'
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']:
                continue
            self.component_call(self.services, port_name, port_dict[port_name], init_mode, 0)

        # -- get plasma state files into driver work directory
        self.services.stage_state()

        # -- post init processing: stage output
        self.services.stage_output_files(0, self.OUTPUT_FILES)

        # -- iteration loop
        nstep = int(self.services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print(("number of interation :", nstep))

        for k in range(nstep):
            print(('\n'+72*"="))
            print(('= model driver: iteration number = %d\n' % k))

            time_id = k

            self.services.update_time_stamp(time_id)

            # call step for each component
            for port_name in port_names:
                if port_name in ['INIT', 'DRIVER'] + post_names:
                    continue
                self.component_call(self.services, port_name, port_dict[port_name], 'step', time_id)

            self.services.stage_state()
            self.services.stage_output_files(time_id, self.OUTPUT_FILES)

        # -- port in post process
        for port_name in post_names:
            self.component_call(self.services, port_name, port_dict[port_name], 'step', time_id)

        self.services.stage_state()
        self.services.stage_output_files(time_id, self.OUTPUT_FILES)

        self.services.update_state()

        # -- call finalize on each component
        print('')
        print((72*"="))
        print('= modeleq driver: finalize')
        print('')
        for port_name in port_names:
            if port_name in ['INIT', 'DRIVER']:
                continue
            self.component_call(self.services, port_name, port_dict[port_name], 'finalize', time_id)

    def finalize(self, timeid=0):
        print('* MODELEQ DRIVER FINALIZE')

        sym_root = self.services.get_config_param('SIM_ROOT')
        outfile = os.path.join(sym_root, 'RESULT')
        print(outfile)
        f = open(outfile, "w")
       #f.write("%d scan_id\n"%0)
        f.write("%d\n" % 0)
        f.close()

        dir_summary = getattr(self, "SUMMARY", "")
        if dir_summary != "":
            if not os.path.exists(dir_summary):
                os.makedirs(dir_summary)
            dir_state = self.services.get_config_param('STATE_WORK_DIR')
            for filename in glob.glob(os.path.join(dir_state, "*.*")):
                print(filename)
                shutil.copy(filename, dir_summary)

    def component_call(self, services, port_name, comp, mode, time):
        comp_mode_string = port_name + ' ' + mode
        t0 = timer.time()
        print(comp_mode_string)
        try:
            services.call(comp, mode, time)
        except Exception:
            services.exception(comp_mode_string + ' failed')
            raise Exception(comp_mode_string + ' failed')
        else:
            print((comp_mode_string + ' finished'))
        t1 = timer.time()
        print(("WALL TIMER: %s %6.3f" % (port_name, t1-t0)))
