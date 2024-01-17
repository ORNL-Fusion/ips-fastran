"""
 -----------------------------------------------------------------------
 timetrace constraint component
 -----------------------------------------------------------------------
"""
import numpy as np
import netCDF4
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate
from fastran.state.timetrace import Timetrace


class constraint_timetrace(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> constraint_timetrace.init() started')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

    def step(self, timeid=0):
        print('>>> constraint_timetrace.step() started')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        next_instate_file = self.services.get_config_param('NEXT_INSTATE')
        next_state_file = self.services.get_config_param('NEXT_STATE')
        next_eqdsk_file = self.services.get_config_param('NEXT_EQDSK')

        # -- load instate and ps
        instate = Instate(cur_instate_file)

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        # -- configuration variables
        tmin = float(self.services.sim_conf['ITERATION_LOOP']['TMIN'])
        tmax = float(self.services.sim_conf['ITERATION_LOOP']['TMAX'])
        dt = float(self.services.sim_conf['ITERATION_LOOP']['DT'])
        f_timetrace = getattr(self, 'TIMETRACE', 'timetrace.nc')
        interpolation_method = getattr(self, 'INTERPOLATION_METHOD', 'linear')
        update_next = getattr(self, 'UPDATE_NEXT', 'disabled')
        exclude = getattr(self, 'EXCLUDE', '')

        # -- time stamp
        timenow = tmin + timeid * dt

        t0 = timenow
        t1 = timenow + dt
        print (f't0, t1 = {t0, t1}')

        instate['t0'] = [t0]
        instate['t1'] = [t1]

        ps['t0'] = t0 * 1.e-3
        ps['t1'] = t1 * 1.e-3

        # -- exclude profiles from timetrace to instate/ps
        if timeid > 0:
            for key in exclude.split():
                instate[f'trace_{key}'] = [0]

        instate.from_timetrace(f_timetrace, timenow, interpolation_method='linear')
        instate.particle_balance()
        instate.scale()

        instate.to_ps_profile(ps)
        if int(getattr(self, 'TRACE_NB', '0')): instate.to_ps_pnbi(ps) 
        # if int(getattr(self, 'TRACE_NB', '0')): instate.to_ps_nb(ps) 
        instate.write(cur_instate_file)
        ps.store(cur_state_file)

        # -- next state
        if update_next == 'enabled':
            instate.from_timetrace(f_timetrace, timenow + dt, interpolation_method='linear')
            instate.particle_balance()
            instate.scale()

            instate.to_ps_profile(ps)
            # if int(getattr(self, 'TRACE_NB', '0')): instate.to_ps_nb(ps)
            instate.write(next_instate_file)
            ps.store(next_state_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('>>> constraint_timetrace.finalize() called')

        # include = getattr(self, "INCLUDE", "").split()
        # if "NB" in include:
        #     pinj = timetrace.slice("pinj", timeid)
        #     print("pinj =", pinj)
        #     instate["power_nbi"] = pinj
        #     ps["power_nbi"] = pinj
