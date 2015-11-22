#!/usr/bin/env python

#------------------------------------------------------------------------------
# data class for fastran : efit
# JM
#------------------------------------------------------------------------------

import os,shutil
from numpy import *
import  netCDF4
from Namelist import Namelist
import zefitutil

FASTRAN_ROOT = os.environ["FASTRAN_ROOT"]
DIR_BIN = os.path.join(FASTRAN_ROOT,"bin")
pstool = os.path.join(DIR_BIN,"pstool")

class efitdata():

    def __init__(self,file_geqdsk,nrho_eq=101,nth_eq=101,nrho_eq_geo=101,f_geq=""):

        if f_geq:
            print 'f_geq',f_geq
            self.data = zefitutil.readg(file_geqdsk) 
            self.nc = netCDF4.Dataset(f_geq,"r",format='NETCDF4')
        else: 
            self.data = zefitutil.readg(file_geqdsk) 
    
            if file_geqdsk != "geqdsk":
                shutil.copyfile(file_geqdsk,"geqdsk")       
    
            inps = Namelist()
            inps["inps"]["nrho_eq"] = [nrho_eq] 
            inps["inps"]["nth_eq"] = [nth_eq] 
            inps["inps"]["nrho_eq_geo"] = [nrho_eq_geo]
            inps.write("inps")
    
            os.system(pstool+" geqdsk 1.0e-6 >& pstool_geqdsk.log")
    
            self.nc = netCDF4.Dataset("ps.nc","r",format='NETCDF4')

    def __getitem__(self, key):

        if key in self.data.keys(): return self.data[key]
        elif key in self.nc.variables: return self.nc.variables[key]
        else: return []

if __name__ == "__main__":

    ps = efitdata("geqdsk_test")
    print ps["psipol"][:]


