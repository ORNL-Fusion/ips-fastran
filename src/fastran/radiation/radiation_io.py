"""
 -----------------------------------------------------------------------
 utils for genray IO
 -----------------------------------------------------------------------
"""
import os
import numpy as np
from scipy.interpolate import RectBivariateSpline
import netCDF4


class radiation_io():
    def __init__(self, dir_data='.'):
        self.dir_data = dir_data
        pass
     
    tiny = np.finfo(np.float64).tiny

    def log10_with_floor(self, x):
        """Return the log of x if x > 0, and otherwise return the log of the smallest representable float."""
        return np.log10(np.maximum(x, self.tiny))
    
    def read(self, impurity_name):
        #-- read genray output
        ncfile = os.path.join(self.dir_data, f'{impurity_name}.nc')
        data = netCDF4.Dataset(ncfile, 'r', format='NETCDF4')
        data_Lz = data.variables['coronal_Lz'][:, :]
        dim_ne_tau = data.variables['dim_ne_tau'][:]
        dim_ne = data.variables['dim_electron_density'][:]
        dim_te = data.variables['dim_electron_temp'][:]
        print(dim_ne_tau)
        print(dim_ne[0], dim_ne[-1])
        print(dim_te[0], dim_te[-1])
        """Data is log scale so must take log off all quantities before interpolation"""
        self.Lz = RectBivariateSpline(self.log10_with_floor(dim_te), self.log10_with_floor(dim_ne), self.log10_with_floor(data_Lz)) 
       
    def __call__(self, te, ne, nte_tau=0.5e17):
        '''Return 10 to the power of result'''
        return np.power(10, self.Lz(self.log10_with_floor(te), self.log10_with_floor(ne), grid=False)) #must be temp in eV, density in m^-3
