#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 driver for fastran scenario modeling
 -----------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *

from component import Component

class fastran_driver(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    # ------------------------------------------------------------------
    # init function

    def init(self, timestamp=0):

        print 80*'*'
        print 'FASTRAN DIRVER INIT'

    # ------------------------------------------------------------------
    # step function

    def step(self, timestamp=0):

        services = self.services

        #-- stage input and plasma state files

        services.stage_input_files(self.INPUT_FILES)
        services.stage_plasma_state()

        #-- get list of ports

        ports = services.getGlobalConfigParameter('PORTS')
        port_names = ports['NAMES'].split()
        print 'PORTS =', port_names
        port_dict = {}

        #-- instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name)
            if(port == None):
                print 'Error accessing %s component'%port_name
                raise
            port_dict[port_name] = port

        #-- simulation mode

        sim_mode = services.getGlobalConfigParameter('SIMULATION_MODE')
        print 'SIMULATION_MODE =', sim_mode

        #-- initial iteration id

        if not hasattr(self,'SUB_WORKFLOW'): self.SUB_WORKFLOW = ''

        if self.SUB_WORKFLOW  == '':
            kstep_prefix = ''
        else:
            kstep_prefix = "iteration_%d_"%timestamp

        t0 = kstep_prefix + "0"

        print '*** t0=',t0, kstep_prefix

        #-- initialize components in PORTS list for startup or restart

        init_mode = 'init'
        if sim_mode == 'RESTART' : init_mode = 'restart'

        if self.SUB_WORKFLOW  == '':
            for port_name in port_names:
                if port_name in ['INIT','DRIVER']: continue
                print port_name, 'init ********'
                self.component_call(services,port_name,port_dict[port_name],init_mode,t0)

        #-- get plasma state files into driver work directory

        services.stage_plasma_state()

        #-- post init processing: stage output

        services.stage_output_files(t0, self.OUTPUT_FILES)

        #-- initial equilibrium

        # if "EQ" in port_names:
        #     print "INITIAL EQUILIRBIUM"
        #     services.call(port_dict["EQ"], "control", {"PRESSURE":"kinetic"})

        #-- steady-state solution procedures, timeloop refers to nonlinear iteration

        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print "number of interation :",nstep

        for kstep in range(nstep):

            # t = kstep_prefix+"%d"%kstep
            t = kstep

            print ''
            print 72*"="
            print '= fastran driver: iteration number = ', kstep
            print ''
            services.update_time_stamp(t)

            # services.stage_plasma_state()
            # self.pre_step_logic(float(t))
            # services.update_plasma_state()

            for port_name in port_names:

                #if port_name in ["INIT","DRIVER"]: continue
                #if port_name in ["EC","IC","NB"]:
                #    if kstep%5 != 0: continue

                if port_name in ['EQ','EC','IC','HC','LH','NB','TR','TR0','BC','MONITOR']:
                    self.component_call(services,port_name,port_dict[port_name],'step',t)

            services.stage_plasma_state()
            services.stage_output_files(t, self.OUTPUT_FILES)

        #-- call finalize on each component

        if 'POST' in port_names:

            nstep_post = int(services.sim_conf["ITERATION_LOOP"]["NSTEP_POST"])
            for kstep in range(nstep,nstep+nstep_post):
                print ''
                print 72*"="
                print '= POST PROCESS: iteration number = ', kstep
                print ''
                # t = kstep_prefix+"%d"%kstep
                t = kstep
                self.component_call(services,'POST',port_dict['POST'],'step',t)
               #self.component_call(services,'EQ',port_dict['EQ'],'step',t)

        if self.SUB_WORKFLOW  == '':

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
        print 'FASTRAN DIRVER FINALIZE'
        sym_root = self.services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
        f.write("%6.3f 0.0\n")
        f.write("%6.3f 0.0\n")
        f.close()

    # ----------------------------------------------------------------------
    # "Private" driver methods

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
