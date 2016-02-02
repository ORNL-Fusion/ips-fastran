#!/usr/bin/env python

#------------------------------------------------------------------------------
# data class for fastran : efit
#------------------------------------------------------------------------------

import os,shutil
from numpy import *
import  netCDF4
from Namelist import Namelist
import zefitutil
import zplasmastate

class efitdata():

    def __init__(self,file_geqdsk,nrho_eq=101,nth_eq=101,nrho_eq_geo=101,f_geq=""):

        if f_geq:
            print 'f_geq',f_geq
            self.data = zefitutil.readg(file_geqdsk) 
            self.ps = zplasmastate.zplasmastate('ips',1)
            self.ps.read(f_geq)

        else: 
            self.data = zefitutil.readg(file_geqdsk) 
            self.ps = zplasmastate.zplasmastate('ips',1)
            self.ps.init_from_geqdsk (file_geqdsk,nrho=129,nth=101)

    def __getitem__(self, key):

        if key in self.data.keys(): return self.data[key]
        else: return self.ps[key]

if __name__ == "__main__":

    ps = efitdata("geqdsk")
    print ps["psipol"][:]


