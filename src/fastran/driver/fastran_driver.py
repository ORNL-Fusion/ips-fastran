"""
 -----------------------------------------------------------------------
 driver for fastran scenario modeling
 -----------------------------------------------------------------------
"""

import os
import glob
import shutil
from ipsframework import Component


class fastran_driver(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print(80*'*')
        print('FASTRAN DIRVER INIT')

    def step(self, timeid=0):
        # -- stage input and plasma state files
        self.services.stage_input_files(self.INPUT_FILES)
        self.services.stage_state()

        # -- get list of ports
        ports = self.services.get_config_param('PORTS')
        port_names = []
        for port_name in ports['NAMES'].split():
            if port_name not in port_names: port_names.append(port_name)

        print('PORTS (all) =', ports['NAMES'].split())
        print('PORTS =', port_names)

        if 'PREPROCESS' in ports:
            PREPROCESS = ports['PREPROCESS'].split()
        else:
            PREPROCESS = []
        if 'POSTPROCESS' in ports:
            POSTPROCESS = ports['POSTPROCESS'].split()
        else:
            POSTPROCESS = []

        # -- instantiate components in port_names list, except DRIVER itself
        port_dict = {}
        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = self.services.get_port(port_name)
            if(port == None):
                raise Exception('Error accessing %s component'%port_name)
            port_dict[port_name] = port

        # -- simulation mode
        sim_mode = self.services.get_config_param('SIMULATION_MODE')
        print('SIMULATION_MODE =', sim_mode)

        #-- initial iteration id
        self.SUB_WORKFLOW = getattr(self, 'SUB_WORKFLOW', '')
        if self.SUB_WORKFLOW  == '':
            kstep_prefix = ''
        else:
            kstep_prefix = "iteration_%d_"%timeid

        t0 = kstep_prefix + "0"
        print('*** t0=', t0, kstep_prefix)

        # -- initialize components in PORTS list for startup or restart
        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        if self.SUB_WORKFLOW  == '':
            for port_name in port_names:
                if port_name in ['INIT', 'DRIVER']: continue
                print(port_name, 'init ********')
                self.services.call(port_dict[port_name], init_mode, t0)

        # -- get plasma state files into driver work directory
        self.services.stage_state()

        # -- post init processing: stage output
        self.services.stage_output_files(t0, self.OUTPUT_FILES, save_plasma_state=False)

        # -- steady-state solution procedures, timeloop refers to nonlinear iteration
        nstep = int(self.services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print("number of interation :", nstep)

        # -- pre process
        for port_name in PREPROCESS:
            self.component_call(port_name, port_dict[port_name], 'step', timeid)

        # -- main iteration
        for kstep in range(timeid, timeid+nstep):
            # t = kstep_prefix+"%d"%kstep
            t = kstep

            print('')
            print(72*"=")
            print('= fastran driver: iteration number = ', kstep)
            print('')
            self.services.update_time_stamp(t)

            for port_name in port_names:
                if port_name in ["INIT", "DRIVER"]: continue
                if port_name in PREPROCESS: continue
                if port_name in POSTPROCESS: continue

                self.component_call(port_name, port_dict[port_name], 'step', t)

            self.services.stage_state()
            self.services.stage_output_files(t, self.OUTPUT_FILES, save_plasma_state=False)

        # -- post process
        if POSTPROCESS:
            nstep_post = int(self.services.sim_conf["ITERATION_LOOP"]["NSTEP_POST"])
            for kstep in range(nstep, nstep+nstep_post):
                # t = kstep_prefix+"%d"%kstep
                t = kstep
                print('')
                print(72*"=")
                print('= POST PROCESS: iteration number = ', t)
                print('')
                self.services.update_time_stamp(t)

                for port_name in POSTPROCESS:
                    self.component_call(port_name, port_dict[port_name], 'step', t)

        # -- call finalize on each component
        if self.SUB_WORKFLOW  == '':
            print('')
            print(72*"=")
            print('= fastran driver: finalize')
            print('')
            for port_name in port_names:
                if port_name in ['INIT', 'DRIVER']: continue
                self.component_call(port_name, port_dict[port_name], 'finalize', t)

    def finalize(self, timeid=0):
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

    def component_call(self, port_name, comp, mode, time):
        comp_mode_string = port_name + ' ' + mode
        print (comp_mode_string)
        try:
            self.services.call(comp, mode, time)
        except Exception:
            self.services.exception(comp_mode_string + ' failed')
            raise Exception("component call")
        else:
            print(comp_mode_string + ' finished')
