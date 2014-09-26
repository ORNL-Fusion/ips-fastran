#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 get ech parameters
 last modified Jun, 2013
 ----------------------------------------------------------------------
"""

import os,sys,pickle
from numpy import *
from scipy import interpolate
import Namelist
import pmds

import netCDF4
import zinterp
import zplot

import xfastran_env

def boxsmooth(y,n):
    ys = zeros(len(y))
    for i in range(n,nlen(y)-n):
        tmp = 0
        for k in range(2*n+1):
            tmp+=y[i-n-k]
        tmp = tmp/(2*n+1.0)
        ys[i] = ys
    return ys

def get_ech(shot,tmin=2000,tmax=6000,dt=20,outdir='.',iplot=False):

    #------------------------------------------------------------------
    # connect mds and open rf tree

    print 'connecting to mds server'
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    print 'open rf tree'
    pmds.mdsopen('RF',shot)

    #------------------------------------------------------------------
    # get data

    print 'get data'
    num_antennas = pmds.mdsvalue('\TOP.ECH:NUM_ANTENNAS')
    num_systems = pmds.mdsvalue('\TOP.ECH:NUM_SYSTEMS')
    ectime = pmds.mdsvalue('\TOP.ECH:ECTIME')
    ectimes = pmds.mdsvalue('\TOP.ECH:ECTIMES')

    times = arange(tmin,tmax,dt)
    num_times = len(times)

    gname = num_systems*[""]

    pinj = zeros((num_systems,num_times))
    aziang = zeros((num_systems,num_times))
    polang = zeros((num_systems,num_times))
    frequency = zeros(num_systems)
    dispersion = zeros(num_systems)
    launch_r = zeros(num_systems)
    launch_z = zeros(num_systems)

    for k in range(num_systems):

        id = k+1
        gname[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.GYROTRON:NAME'%id)
        frequency[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.GYROTRON:FREQUENCY'%id)
        dispersion[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:DISPERSION'%id)
        launch_r[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:LAUNCH_R'%id)
        launch_z[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:LAUNCH_Z'%id)
        aziang_a = pmds.mdsvalue('\EC%sAZIANG'%gname[k][0:3])
        polang_a = pmds.mdsvalue('\EC%sPOLANG'%gname[k][0:3])
        pinj_a = pmds.mdsvalue('\EC%sFPWRC'%gname[k][0:3])

        print 'processing: %d, %s'%(k,gname[k])

        for i in range(num_times):
            ind = (ectime >= times[i]-0.5*dt ) & (ectime <= times[i]+0.5*dt )
            pinj[k,i] = average(pinj_a[ind])
            if pinj[k,i] <= 3.5e4: pinj[k,i]=0.0

        aziang_spl = zinterp.zinterp(ectimes,aziang_a)
        polang_spl = zinterp.zinterp(ectimes,polang_a)

        aziang[k] = aziang_spl(times)
        polang[k] = polang_spl(times)

    #------------------------------------------------------------------
    # disconnect mds

    pmds.mdsdisconnect()

    #------------------------------------------------------------------
    # check plot

    if iplot:
        p = zplot.zplot("ech_%06d.pdf"%shot,size_x=8.0,size_y=4.0)
        for k in range(num_systems):
            zplot.plot_s(111,1.0e-3*times,1.0e-6*pinj[k,:]
                 ,[0.0,7.0,1.0],[0.0,1.0,0.2]
                ,iline='line',isym='',xlab='time',ylab='pinj')
            p.savefig()
        p.close()

    #------------------------------------------------------------------
    # output

    ech = { 
        "num_systems" : num_systems,  
        "gname"       : gname,
        "num_times"   : num_times,  
        "times"       : times,      # [msec]
        "frequency"   : frequency,  # [Hz]
        "dispersion"  : dispersion, # []
        "launch_r"    : launch_r,   # [m]
        "launch_z"    : launch_z,   # [m]
        "aziang"      : aziang,     # [degree]
        "polang"      : polang,     # [degree]
        "pinj"        : pinj        # [w]
    }

    #------------------------------------------------------------------
    # dump to netcdf

    print 'dump to '+outdir+'/ech_%06d.nc'%shot

    nc = netCDF4.Dataset(outdir+'/ech_%06d.nc'%shot,'w'
             ,format=xfastran_env.netcdf_form)
    
    dim_id   = nc.createDimension('id'  , ech["num_systems"])
    dim_time = nc.createDimension('time', ech["num_times"  ])
    
    id   = nc.createVariable('id','i4',('id',))
    time = nc.createVariable('time','f8',('time',))

    frequency  = nc.createVariable('frequency','f8',('id',))  
    dispersion = nc.createVariable('dispersion','f8',('id',)) 
    launch_r   = nc.createVariable('launch_r','f8',('id',))      
    launch_z   = nc.createVariable('launch_z','f8',('id',))      
    aziang     = nc.createVariable('aziang','f8',('id','time',))      
    polang     = nc.createVariable('polang','f8',('id','time',))       
    pinj       = nc.createVariable('pinj','f8',('id','time',))        
    
    id[:]   = range(ech["num_systems"])
    time[:] = ech["times"]

    frequency[:]  = ech["frequency"] 
    dispersion[:] = ech["dispersion"]
    launch_r[:]   = ech["launch_r"]  
    launch_z[:]   = ech["launch_z"]  
    aziang[:,:]   = ech["aziang"]         
    polang[:,:]   = ech["polang"]          
    pinj[:,:]     = ech["pinj"]           
    
    nc.close()

    #------------------------------------------------------------------
    # return

    return ech

if __name__ == "__main__":

    shot=153646
    ech = get_ech(shot)




