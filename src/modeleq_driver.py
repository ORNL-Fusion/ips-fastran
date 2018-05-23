#! /usr/bin/env python

"""
 fastran model equilibrium driver
"""

import os
import time as timer

from component import Component
from Namelist import Namelist

class modeleq_driver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0):

        print '* DIRVER INIT'

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

        post_names = ports['POSTS'].split()
        print 'POSTS =', post_names

        #-- instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name)
            if(port == None):
                raise Exception('Error accessing %s component'%port_name)
            port_dict[port_name] = port

        #-- initialize components in PORTS list for startup

        init_mode = 'init'
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue
            self.component_call(services,port_name,port_dict[port_name],init_mode,0)

        #-- get plasma state files into driver work directory

        services.stage_plasma_state()

        #-- post init processing: stage output

        services.stage_output_files(0, self.OUTPUT_FILES)

        #-- iteration loop

        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print "number of interation :",nstep

        for k in range(nstep):

            print '\n'+72*"="
            print '= model driver: iteration number = %d\n'%k

            time_id = k

            services.update_time_stamp(time_id)

            # call step for each component

            for port_name in port_names:
                if port_name in ['INIT','DRIVER'] + post_names: continue
                self.component_call(services,port_name,port_dict[port_name],'step',time_id)

            services.stage_plasma_state()
            services.stage_output_files(time_id, self.OUTPUT_FILES)

        #-- port in post process

        for port_name in post_names:
            self.component_call(services,port_name,port_dict[port_name],'step',time_id)

        services.stage_plasma_state()
        services.stage_output_files(time_id, self.OUTPUT_FILES)

        services.update_plasma_state()

        #-- call finalize on each component

        print ''
        print 72*"="
        print '= modeleq driver: finalize'
        print ''
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue
            self.component_call(services, port_name, port_dict[port_name], 'finalize', time_id)

    def finalize(self, timestamp = 0):

        print 'DRIVER FINALIZE'

        services = self.services
        sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        print outfile
        f = open(outfile,"w")
        f.write("%d scan_id\n"%0)
        f.close()
        pass

    def component_call(self, services, port_name, comp, mode, time):

        comp_mode_string = port_name + ' ' + mode
        t0 = timer.time()
        print (comp_mode_string)
        try:
            services.call(comp, mode, time)
        except Exception:
            services.exception(comp_mode_string + ' failed')
            raise
        else:
            print comp_mode_string + ' finished'
        t1 = timer.time()
        print "WALL TIMER: %s %6.3f"%(port_name, t1-t0)
        return