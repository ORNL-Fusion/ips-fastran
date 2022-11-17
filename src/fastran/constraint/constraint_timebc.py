"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""
import numpy as np
import netCDF4
from scipy.interpolate import interp1d
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate
from fastran.constraint.constraint_timebc_io import Timetrace


class constraint_timebc(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('constraint_timebc.init() started')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- timebc data
        timebc = Namelist('intimebc')
        self.data = {}

        for key in timebc['timebc'].keys():
            print(key)
            if key.split('_')[0] == 'TIME':
               vname = key.split('TIME_')[1]
               self.data[vname] = {'t': timebc['timebc'][key], 'y': timebc['timebc'][vname]}
        print(self.data)

    def _slice(self, key, time):
        return interp1d(self.data[key]['t'], self.data[key]['y'], axis=0, fill_value=(self.data[key]['y'][0], self.data[key]['y'][-1],), bounds_error=False)(time)

    def step(self, timeid=0):
        print('>>> constraint_timebc.step() started')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- from timetrace
        tmin = float(self.services.sim_conf["ITERATION_LOOP"]["TMIN"])
        dt = float(self.services.sim_conf["ITERATION_LOOP"]["DT"])
        timenow = tmin + timeid*dt
        print('timenow = ', timenow)
        include = getattr(self, "INCLUDE", "").split()

        instate = Instate(cur_instate_file)

        for key in self.data.keys():
            y = self._slice(key, timenow)
            print('TIME BC: {} = {}'.format(key, y))
            instate[key] = [y]

        instate.write(cur_instate_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('constraint_timebc.finalize() called')

