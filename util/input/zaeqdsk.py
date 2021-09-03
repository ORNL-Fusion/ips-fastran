#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 collect efit geqdsks 
 last modified Jun, 2013
 ----------------------------------------------------------------------
"""

import os
import shutil
import pickle
import glob
from numpy import *

import netCDF4
import efit_eqdsk

netcdf_form = 'NETCDF3_CLASSIC'

def get_shot_time(fname):

    shot = int(os.path.basename(fname).split('.')[0][-6:])
    time = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    return (shot,time)

def get_efit(shot,tmin=2000,tmax=6000,outdir='.',efitdir='.'):

    afiles = glob.glob(os.path.join(efitdir,'a%06d.?????'%shot))
    afiles.sort()

    times = []
    for afile in afiles:
        time = get_shot_time(afile)[-1]
        if time<tmin or time>tmax: continue
        times.append(time)

    ntimes = len(times)

    key_list = [
        "qaxis","rout","r0","a0","ip","b0","betat","betan",
        "rmajor","kappa","delta_u","delta_l","betap","li","peak","q95","qmin","rho_qmin"
    ]

    data = {}
    for key in key_list: data[key] = zeros(ntimes)

    for k in range(ntimes):

        afile = os.path.join(efitdir,'a%06d.%05d'%(shot,times[k]))
        print 'processing',afile

        a  = efit_eqdsk.reada(afile)
   
        for key in key_list:
            data[key][k] = a[key]

    print 'dump to '+outdir+'/aeq_%06d.nc'%shot

    nc = netCDF4.Dataset(outdir+'/aeq_%06d.nc'%shot,'w' ,format=netcdf_form)

    
    dim_time = nc.createDimension('time', ntimes)
    
    var_time   = nc.createVariable('time','f8',('time',))
    var_time[:] = times

    data_nc = {}
    for key in key_list: 
        data_nc[key] = nc.createVariable(key,'f8',('time',))

    for key in key_list:
        data_nc[key][:] = data[key] 

    nc.close()


if __name__ == "__main__":

    shot= 153646
    efitdir = '/u/parkjm/013/20132204/153646/kin02'
    tmp = get_efit(shot,3700,3800,outdir='.',efitdir=efitdir)




