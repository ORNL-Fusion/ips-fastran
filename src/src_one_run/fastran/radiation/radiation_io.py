"""
 -----------------------------------------------------------------------
 utils for genray IO
 -----------------------------------------------------------------------
"""
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import netCDF4


class radiation_io():
    def __init__(self, dir_data='.'):
        self.dir_data = dir_data
        pass
    
    def read(self, impurity_name):
        #-- read genray output
        ncfile = os.path.join(self.dir_data, f'{impurity_name}.nc')
        data = netCDF4.Dataset(ncfile, 'r', format='NETCDF4')
        data_Lz = data.variables['equilibrium_Lz'][:, :, :]
        dim_ne_tau = data.variables['dim_ne_tau'][:]
        dim_ne = data.variables['dim_electron_density'][:]
        dim_te = data.variables['dim_electron_temp'][:]
        print(dim_ne_tau)
        print(dim_ne[0], dim_ne[-1])
        print(dim_te[0], dim_te[-1])
        self.Lz = RegularGridInterpolator((dim_te, dim_ne, dim_ne_tau), data_Lz) 
    
    def __call__(self, te, ne, nte_tau=0.5e17):
        return self.Lz((te, ne, nte_tau))