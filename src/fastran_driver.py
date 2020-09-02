"""
 -----------------------------------------------------------------------
 driver for fastran scenario modeling
 -----------------------------------------------------------------------
"""

import os
from component import Component

class fastran_driver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print(80*'*')
        print('FASTRAN DIRVER INIT')

    def step(self, timeid=0):
        #--- entry
        services = self.services

        #-- stage input and plasma state files
        services.stage_input_files(self.INPUT_FILES)
        services.stage_plasma_state()

        #-- get list of ports
        ports = services.getGlobalConfigParameter('PORTS')
        #port_names = ports['NAMES'].split()
        port_names = []
        for port_name in ports['NAMES'].split():
            if port_name not in port_names: port_names.append(port_name)

        print('PORTS (all) =', ports['NAMES'].split())
        print('PORTS =', port_names)

        if 'PREPROCESS' in ports:
            PREPROCESS = ports['PREPROCESS'].split()
        else:
            PREPROCESS = []

        #-- instantiate components in port_names list, except DRIVER itself
        port_dict = {}
        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name)
            if(port == None):
                raise Exception('Error accessing %s component'%port_name)
            port_dict[port_name] = port

        #-- simulation mode
        sim_mode = services.getGlobalConfigParameter('SIMULATION_MODE')
        print('SIMULATION_MODE =', sim_mode)

        #-- initial iteration id
        self.SUB_WORKFLOW = getattr(self, 'SUB_WORKFLOW', '')
        if self.SUB_WORKFLOW  == '':
            kstep_prefix = ''
        else:
            kstep_prefix = "iteration_%d_"%timeid

        t0 = kstep_prefix + "0"

        print('*** t0=',t0, kstep_prefix)

        #-- initialize components in PORTS list for startup or restart
        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        if self.SUB_WORKFLOW  == '':
            for port_name in port_names:
                if port_name in ['INIT', 'DRIVER']: continue
                print(port_name, 'init ********')
                services.call(port_dict[port_name], init_mode, t0)
                #self.component_call(services, port_name, port_dict[port_name], init_mode, t0)

        #-- get plasma state files into driver work directory
        services.stage_plasma_state()

        #-- post init processing: stage output
        services.stage_output_files(t0, self.OUTPUT_FILES)

        #-- steady-state solution procedures, timeloop refers to nonlinear iteration
        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print("number of interation :", nstep)

        #-- pre process
        for port_name in PREPROCESS:
            self.component_call(services, port_name, port_dict[port_name], 'step', 0)

        #-- main iteration
        for kstep in range(timeid, timeid+nstep):
            # t = kstep_prefix+"%d"%kstep
            t = kstep

            print('')
            print(72*"=")
            print('= fastran driver: iteration number = ', kstep)
            print('')
            services.update_time_stamp(t)

            for port_name in port_names:

                if port_name in PREPROCESS: continue
                if port_name in ['EQ', 'EC', 'IC', 'HC', 'LH', 'NB', 'TR', 'BC', 'FEEDBACK', 'MONITOR', 'PEDESTAL']:
                    self.component_call(services, port_name, port_dict[port_name], 'step', t)

            services.stage_plasma_state()
            services.stage_output_files(t, self.OUTPUT_FILES)

        #-- post process
        if 'POST' in port_names:
            nstep_post = int(services.sim_conf["ITERATION_LOOP"]["NSTEP_POST"])
            for kstep in range(nstep, nstep+nstep_post):
                print('')
                print(72*"=")
                print('= POST PROCESS: iteration number = ', kstep)
                print('')
                # t = kstep_prefix+"%d"%kstep
                t = kstep
                self.component_call(services, 'POST', port_dict['POST'], 'step', t)

        #-- call finalize on each component
        if self.SUB_WORKFLOW  == '':
            print('')
            print(72*"=")
            print('= fastran driver: finalize')
            print('')
            for port_name in port_names:
                if port_name in ['INIT', 'DRIVER']: continue
                self.component_call(services, port_name, port_dict[port_name], 'finalize', t)

    def finalize(self, timeid = 0):
        print(80*'*')
        print('FASTRAN DIRVER FINALIZE')
        sym_root = self.services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f\n"%0.0)
        f.write("%6.3f\n"%0.0)
        f.close()

    def component_call(self, services, port_name, comp, mode, time):
        comp_mode_string = port_name + ' ' + mode
        print (comp_mode_string)
        try:
            services.call(comp, mode, time)
        except Exception:
            services.exception(comp_mode_string + ' failed')
            raise
        else:
            print(comp_mode_string + ' finished')
