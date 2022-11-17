from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate

import numpy as np
import netCDF4
from Namelist import Namelist

from scipy.interpolate import interp1d

class Timetrace():
    def __init__(self, fname, start_time=-1, include=[]):
        timetrace = Namelist(fname)
        self.data = {}

        for key in timetrace['timetrace'].keys():
            print(key)
            if key.split('_')[0] == 'TIME':
               vname = key.split('TIME_')[1]
               self.data[vname] = {'t': timetrace['timetrace'][key], 'y': timetrace['timetrace'][vname]}
        print(self.data)

    def _slice(self, key, time):
        return interp1d(self.data[key]['t'], self.data[key]['y'], axis=0, fill_value=(self.data[key]['y'][0], self.data[key]['y'][-1],), bounds_error=False)(time)

    def slice(self, key, i):
        return self.get(key, i)

if __name__=="__main__":
    timetrace = Timetrace("intimetrace")
    print( timetrace._slice('IP', 0.3))
    for t in np.linspace(0.0, 0.5, 100):
        print(t, timetrace._slice('IP', t))
