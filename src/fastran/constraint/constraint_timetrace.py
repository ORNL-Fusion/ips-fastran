"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate
from fastran.state.timetrace import Timetrace

import numpy as np
import netCDF4
from Namelist import Namelist


class constraint_timetrace(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('constraint_timetrace.init() started')

        #--- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

    def step(self, timeid=0):
        print('constraint_timetrace.step() started')

        #--- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #--- stage plasma state files
        self.services.stage_state()

        #--- get plasma state file name
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        next_state_file = self.services.get_config_param('NEXT_STATE')
        next_eqdsk_file = self.services.get_config_param('NEXT_EQDSK')

        #-- from timetrace
        tmin = int(self.services.sim_conf["TIME_LOOP"]["TMIN"])
        tmax = int(self.services.sim_conf["TIME_LOOP"]["TMAX"])

        f_timetrace = self.TIMETRACE
        timetrace = Timetrace(f_timetrace, start_time=tmin, end_time=tmax)

        ip = timetrace.get("ip", timeid)
        b0 = timetrace.get("bt", timeid)
        r0 = timetrace.get("r0", timeid)
        print('ip = {}, bt={}, r0={}'.format(ip*1.e-6, b0, r0))

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        instate = Instate(cur_instate_file)

        instate["ip"] = [ip*1.e-6]
        instate["b0"] = [b0]
        instate["r0"] = [r0]
        instate["ne"] = timetrace.slice("ne", timeid)
        instate["te"] = timetrace.slice("te", timeid)
        instate["ti"] = timetrace.slice("ti", timeid)
        instate["omega"] = timetrace.slice("omega", timeid)
        instate["zeff"] = 101*[1.6]

        nbdry = timetrace.get("nbdry", timeid)
        instate["nbdry"] = [nbdry]
        instate["rbdry"] = timetrace.get("rbdry", timeid)[:nbdry]
        instate["zbdry"] = timetrace.get("zbdry", timeid)[:nbdry]

        instate["p_eq"] = timetrace.get("p_eq", timeid)
        instate["j_tot"] = timetrace.get("jpar", timeid)*1.e-6

        t0 = timetrace.get_time(timeid)
        t1 = timetrace.get_time(timeid+1)
        print ("t0, t1 =", t0, t1)

        instate["t0"] = [t0]
        instate["t1"] = [t1]

        ps["t0"] = t0*1.e-3
        ps["t1"] = t1*1.e-3

        include = getattr(self, "INCLUDE", "").split()
        if "NB" in include:
            pinj = timetrace.slice("pinj", timeid)
            print("pinj =", pinj)
            instate["power_nbi"] = pinj
            ps["power_nbi"] = pinj

        instate.to_ps(ps)

        instate.write(cur_instate_file)
        ps.store(cur_state_file)

        #--- update plasma state files
        self.services.update_state()

        #--- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('constraint_timetrace.finalize() called')
