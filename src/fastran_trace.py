#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 driver for fastran scenario modeling
 -----------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *

from component import Component

class fastran_trace(Component):

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

        services = self.services

        #-- stage input and plasma state files

        services.stage_input_files(self.INPUT_FILES)
        services.stage_plasma_state()

        #-- get Portal RUNID and save to a file

        # run_id = services.get_config_param('PORTAL_RUNID')
        # sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        # path = os.path.join(sym_root, 'PORTAL_RUNID')
        # runid_file = open(path, 'a')
        # runid_file.writelines(run_id + '\n')
        # runid_file.close()

        #-- get list of ports

        ports = services.getGlobalConfigParameter('PORTS')
        port_names = ports['NAMES'].split()
        print 'PORTS =', port_names
        port_dict = {}
        port_id_list = []

        #-- instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name) 
            if(port == None):
                print 'Error accessing %s component'%port_name
                raise
            port_dict[port_name] = port
            port_id_list.append(port)

        #-- startup or restart?

        sim_mode = services.getGlobalConfigParameter('SIMULATION_MODE')
        print 'SIMULATION_MODE =', sim_mode

        #-- get timeloop for simulation, not used, t = 0

        t = 0

        #-- initialize components in PORTS list for startup or restart

        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services,port_name,port_dict[port_name],init_mode,t)
        
        #-- get plasma state files into driver work directory 

        services.stage_plasma_state()

        #-- post init processing: stage output

        services.stage_output_files(t, self.OUTPUT_FILES)

        #-- initial equilibrium

        # if "EQ" in port_names:
        #     print "INITIAL EQUILIRBIUM"        
        #     services.call(port_dict["EQ"], "control", {"PRESSURE":"kinetic"})

        #-- steady-state solution procedures, timeloop refers to nonlinear iteration

        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print "number of interation :",nstep

        for k in range(nstep):

            iStamp = k 

            print ''
            print 72*"="
            print '= fastran driver: iteration number = ', k
            print ''
            services.update_time_stamp(iStamp)

            # services.stage_plasma_state()
            # self.pre_step_logic(float(t))
            # services.update_plasma_state()

            for port_name in port_names:
                if port_name == 'EQ':
                    self.component_call(services,'EQ',port_dict['EQ'],'step',iStamp)
                elif port_name == 'EC':
                    self.component_call(services,'EC',port_dict['EC'],'step',iStamp)
                elif port_name == 'IC':
                    self.component_call(services,'IC',port_dict['IC'],'step',iStamp)
                elif port_name == 'NB':
                    self.component_call(services,'NB',port_dict['NB'],'step',iStamp)
                elif port_name == 'TR':
                    self.component_call(services,'TR',port_dict['TR'],'step',iStamp)
                elif port_name == 'PED':
                    self.component_call(services,'PED',port_dict['PED'],'step',iStamp)
                elif port_name == 'BC':
                    self.component_call(services,'BC',port_dict['BC'],'step',iStamp)
                elif port_name == 'MONITOR':
                    self.component_call(services,'MOINTOR',port_dict['MONITOR'],'step',iStamp)
                
            services.stage_plasma_state()
            services.stage_output_files(t, self.OUTPUT_FILES)

        #-- call finalize on each component

        print ''
        print 72*"="
        print '= fastran driver: finalize'
        print ''
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services, port_name, port_dict[port_name], 'finalize', t)

    # ------------------------------------------------------------------
    # checkpoint function

    def checkpoint(self, timestamp=0.0):

        print 'fastran_driver.checkpoint() called'
        
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
    # "Private" driver methods

    # Component call - wraps the exception handling for all component calls
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
    
