#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 plasma state, pstool backend
 -----------------------------------------------------------------------
"""
import os
import shutil
import subprocess
from numpy import *
import netCDF4

from Namelist import Namelist

from plasmastate_base import plasmastate_base

class plasmastate(plasmastate_base):

    def __init__(self,lable,level):

        self.pstool = os.environ["PSTOOL"]
        pass

    def init(self, f_ps, **keyargs):

        nml = Namelist()
        for key, value in keyargs.iteritems():
            nml["inpstool"][key] = value 
        nml.write("inpstool")
        print 'start',f_ps,self.pstool
        retcode = subprocess.call([self.pstool, "init", f_ps])
        print 'end'
        self.read(f_ps) 

    def read(self, f_ps, mode='r+'):

        self.data = netCDF4.Dataset(f_ps,mode,format='NETCDF4')
        self.f_ps = f_ps
        self.mode = mode

    def info(self):

        pass

    def close(self):

        self.data.close()
      # retcode = subprocess.call([self.pstool, "rehash", self.f_ps])

    def store(self,f_ps):

        self.close()
        if f_ps != self.f_ps: shutil.copyfile(self.f_ps,f_ps)

    def load_innubeam (self, f_innubeam="innubeam"):

        self.data.close()
        retcode = subprocess.call([self.pstool, "load_nbi", self.f_ps]) #, f_innubeam])
        self.read(self.f_ps) 

    def load_geqdsk (self,fn_geqdsk,bdy_crat=1.e-6):

        #kcur_option = 1
        #rho_curbrk = 0.9 
        self.data.close()
        retcode = subprocess.call([self.pstool, "load_geqdsk", self.f_ps, fn_geqdsk,"%e"%bdy_crat])
        self.read(self.f_ps) 
        
    def init_from_geqdsk (self,fn_geqdsk,nrho=101,nth=101,shot=0,time=0,bdy_crat=1.e-6):

        self.f_ps = "geqdsk.nc"
        nml = Namelist()
        nml["inpstool"]["nrho_eq"] = [nrho] 
        nml["inpstool"]["nth_eq"] = [nth] 
        nml["inpstool"]["nrho_eq_geo"] = [nrho] 
        nml.write("inpstool")
        print self.f_ps, fn_geqdsk
        retcode = subprocess.call([self.pstool, "geqdsk", self.f_ps, fn_geqdsk, "%e"%bdy_crat])
        self.read(self.f_ps) 

    def __getitem__(self, key):

        if key in self.data.variables: return self.data.variables[key]
        else: return []


if __name__ == "__main__":

    ps = plasmastate("ips",1)

    ps.init(
        "ps.nc", 
        global_label = ['ips'],
        nspec_th = [1],
        nspec_beam = [1], 
        nspec_fusion = [0], 
        nspec_rfmin = [0],
        nspec_gas = [1],
        z_ion = [1],
        a_ion = [2],
        z_imp = [6],
        a_imp = [12],
        z_beam = [1],
        a_beam = [2],
        z_fusion = [2],
        a_fusion = [4],
        z_rfmin = [1],
        a_rfmin = [1],
        z_gas = [1],
        a_gas = [2],
        nrho = [101],
    ) 

    ps.load_geqdsk("g147634.04325")

    rho = linspace(0,1.0,100)
    print rho

    ps["Ts"][:] = 3.0
    ps["Ts"][0,:] = rho
    ps["Ti"] = rho
    rho[0] = 100.0
    ps.store("ps.nc")
   
