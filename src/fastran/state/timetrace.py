from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate

import numpy as np
import netCDF4
from Namelist import Namelist

from scipy.interpolate import interp1d

class Timetrace():
    def __init__(self, fname, start_time=-1, end_time=1000000):
        timetrace = netCDF4.Dataset(fname, 'r', format='NETCDF4')
        self.data = {}

        for key in ["time_efit", "time_ne", "time_te", "time_ti", "time_omega", "time_nz"]:
            self.data[key] =  timetrace.variables[key][:]
        for key in ["rho_efit", "rho"]:
            self.data[key] =  timetrace.variables[key][:]
        for key in ["nlim", "nbdry"]:
            self.data[key] =  timetrace.variables[key][:]
        for key in ["ip", "bt", "r0"]:
            self.data[key] =  timetrace.variables[key][:]

        for key in ["rbdry", "zbdry"]:
            self.data[key] =  timetrace.variables[key][:, :]
        for key in ["rlim", "zlim"]:
            self.data[key] =  timetrace.variables[key][:, :]
        for key in ["p_eq", "jpar"]:
            self.data[key] =  timetrace.variables[key][:, :]
        for key in ["ne", "te", "ti", "omega", "nz"]:
            self.data[key] =  timetrace.variables[key][:, :]

        self.data["time_pinj"] =  timetrace.variables["time_nbi"][:]
        for key in ["pinj"]:
            self.data[key] =  timetrace.variables[key][:, :].transpose()

        timetrace.close()

        self.start_time = start_time
        self.end_time = end_time

        self.times = self.data["time_efit"]

        self.istart = np.where(self.times >= self.start_time)[0][0]
        self.iend = np.where(self.times <= self.end_time)[0][-1]

        print("istart = {}, iend = {}".format(self.istart, self.iend))
        print("tstart = {}, tend = {}".format(self.times[self.istart], self.times[self.iend]))

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, key):
        return True if key in self.data else False 

    def get(self, key, i):
        return self.data[key][self.istart+i]

    def _slice(self, key, time):
        return interp1d(self["time_{}".format(key)], self[key], axis=0, fill_value=(self[key][0], self[key][-1],), bounds_error=False)(time)

    def slice(self, key, i):
        return self._slice(key, self.times[self.istart+i])

    def get_time(self, i):
        return self.times[self.istart+i]

if __name__=="__main__":
    timetrace = Timetrace("t153648.nc")
    ne = timetrace["ne"]
    time_ne = timetrace["time_ne"]
    print(time_ne)
    print(ne[0])
    print(ne[1])
    print(timetrace.slice("ne", 2405))
    print(timetrace.slice("ne", 0))
