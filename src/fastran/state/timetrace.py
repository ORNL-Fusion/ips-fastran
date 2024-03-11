import sys
import numpy as np
import netCDF4
from scipy.interpolate import interp1d, interp2d
from Namelist import Namelist


class Timetrace():
    def __init__(self, fname, kind='linear'):
        nc = netCDF4.Dataset(fname, 'r', format='NETCDF4')
        self.data = {}
        self.interpolation = {}

        # equilibrium
        for key in ['time_eq', 'rho_eq']:
            self.data[key] =  nc.variables[key][:]

        for key in ['p_eq', 'j_tot']:
            t = self.data['time_eq']
            x = self.data['rho_eq']
            f = nc.variables[key][:, :]
            self.data[key] = f
            self.interpolation[key] = interp2d(t, x, f.transpose(), kind=kind)

        for key in ['ip', 'bt', 'r0']:
            t = self.data['time_eq']
            f = nc.variables[key][:]
            self.data[key] =  f
            self.interpolation[key] = interp1d(t, f, kind=kind, fill_value=(f[0], f[-1],), bounds_error=False)

        # profile
        var_profile = ['ne', 'te', 'ti', 'omega', 'nz']
        var_nb = ['pe_nb', 'pi_nb', 'density_beam', 'wbeam']
        var_ext = ['p_rad', 'p_ohm']
        for key in var_profile + var_nb + var_ext:
            if key in nc.variables.keys():
                t = nc.variables[f'time_{key}'][:]
                x = nc.variables[f'rho_{key}'][:]
                f = nc.variables[key][:]
                self.data[f'time_{key}'] = t
                self.data[f'rho_{key}'] = x
                self.data[key] = f
                self.interpolation[key] = interp2d(t, x, f.transpose(), kind=kind)

        # limiter
        self.data['rlim'] =  nc.variables['rlim'][:]
        self.data['zlim'] =  nc.variables['zlim'][:]
        self.data['nlim'] =  len(self.data['rlim'])

        # plasma shape
        self.data['nbdry'] =  nc.variables['nbdry'][:]
        self.data['rbdry'] =  nc.variables['rbdry'][:, :]
        self.data['zbdry'] =  nc.variables['zbdry'][:, :]

        # nbi power
        var_nbi = ['pnbi'] #, ['pnbi', 'einj', 'ffull', 'fhalf']
        for key in var_nbi:
            if key in list(nc.variables.keys()):
                self.data[key] = nc.variables[key][:, :]
                self.interpolation[key] = len(self.data[key]) * [None]
                for k in range(len(self.data[key])):
                    t = nc.variables['time_nbi'][:]
                    f = nc.variables[key][k, :]
                    self.interpolation[key][k] = interp1d(t, f, kind=kind, fill_value=(f[0], f[-1],), bounds_error=False) 

        nc.close()

        self.times = self.data['time_eq']
        self.tmin = self.times[0]
        self.tmax = self.times[-1]

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, key):
        return True if key in self.data else False

    def get(self, key, t):
        return self.interpolation[key](t)

    def slice(self, key, t, x):
        return self.interpolation[key](t, x).flatten()

    def slice_list(self, key, t):
        return [self.interpolation[key][k](t) for k in range(len(self.interpolation[key]))]

    def nearest(self, key, t):
        delta = np.abs(self.times - t)
        i = np.argmin(delta)
        return self.data[key][i]
        

if __name__=="__main__":
    # test
    ncfile = sys.argv[1] # convention: t<shot_number>.nc
    timetrace = Timetrace(ncfile)
    ne = timetrace['ne']

    t = 0.5 * (timetrace.tmin + timetrace.tmax)
    rho = np.linspace(0., 1., 20)
    print(timetrace.slice('ne', t, rho) )

    t = np.linspace(timetrace.tmin, timetrace.tmax, 100)
    rho = 0.5
    print(timetrace.slice('ne', t, rho) )
    print(timetrace.get('ip', t))
