#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 get nbi parameters
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

def get_nbi(shot,tmin=0,tmax=10000,dt=20,outdir='.',iplot=False):

    nb_list = ['30L','30R','15L','15R','21L','21R','33L','33R']

    #------------------------------------------------------------------
    # connect mds and open nb tree

    print 'connecting to mds server'
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    print 'open nb tree'
    pmds.mdsopen('NB',shot)

    #------------------------------------------------------------------
    # get data

    print 'get data'

    time_mds = pmds.mdsvalue('DIM_OF(\\PINJ)')
    pinj_mds = {}
    einj_mds = {}
    for src in nb_list:
        pinj_mds[src] = pmds.mdsvalue('\\PINJ_'+src)
        einj_mds[src] = pmds.mdsvalue('.NB'+src+':NBVAC_SCALAR')

    try:
        btilt_15 = pmds.mdsvalue('.NB15L.OANB:BLPTCH_CAD')
    except:
        try:
           btilt_15 = pmds.mdsvalue('.NB15R.OANB:BLPTCH_CAD')
        except:
           btilt_15 = 0.0

    if btilt_15 > 0.0:
        try:
            stilt_15L = pmds.mdsvalue('.NB15L.OANB:SRCPTCH')/60.0
        except:
            stilt_15L = 0.0
        try:
            stilt_15R = pmds.mdsvalue('.NB15R.OANB:SRCPTCH')/60.0
        except:
            stilt_15R = 0.0
    else:
        stilt_15L = 0.0
        stilt_15R = 0.0

    print 'shot :', shot
    print '150  beam   tilt [deg] :', btilt_15
    print '150L source tilt [min] :', stilt_15L
    print '150R source tilt [min] :', stilt_15R

    #------------------------------------------------------------------
    # process data

    times = arange(tmin,tmax,dt)
    num_times = len(times)

    pinj = zeros((len(nb_list),num_times))
    einj = zeros(len(nb_list))
    nbsrc = len(nb_list)*[""]

    for k, src in enumerate(nb_list):

        for i in range(num_times):
            ind = (time_mds >= times[i]-0.5*dt ) & (time_mds <= times[i]+0.5*dt )
            pinj[k,i] = average(pinj_mds[src][ind])
            if pinj[k,i] <= 3.5e4: pinj[k,i]=0.0

        einj[k] = einj_mds[src] 

        nbsrc[k] = src

    #------------------------------------------------------------------
    # disconnect mds

    pmds.mdsdisconnect()

    #------------------------------------------------------------------
    # check plot

    if iplot:
        p = zplot.zplot("nbi_%06d.pdf"%shot,size_x=8.0,size_y=4.0)
        for k, src in enumerate(nb_list):
            zplot.plot_s(111,1.0e-3*times,1.0e-6*pinj[k,:]
                 ,[0.0,7.0,1.0],[0.0,3.0,1.0]
                ,iline='line',isym='',xlab='time',ylab='pinj_%s'%src)
            p.savefig()
        p.close()

    #------------------------------------------------------------------
    # output

    nbi = { 
        "num_srcs"    : len(nb_list),  
        "nbsrc"       : nbsrc,
        "num_times"   : num_times,  
        "times"       : times,      # [msec]
        "pinj"        : pinj,       # [W]
        "einj"        : einj,       # [eV]
        "tilt_15"     : [btilt_15,stilt_15L,stilt_15R] #[deg]
    }

    #------------------------------------------------------------------
    # dump to netcdf

    print 'dump to '+outdir+'/nbi_%06d.nc'%shot

    nc = netCDF4.Dataset(outdir+'/nbi_%06d.nc'%shot,'w'
             ,format=xfastran_env.netcdf_form)
    
    dim_id   = nc.createDimension('id'  , nbi["num_srcs"])
    dim_time = nc.createDimension('time', nbi["num_times"])
    dim_one  = nc.createDimension('one' , 1)
    dim_char = nc.createDimension('char', 3)
    
    #id   = nc.createVariable('id',str,('id',))

    id   = nc.createVariable('id','c',('id','char'))
    time = nc.createVariable('time','f8',('time',))

    einj = nc.createVariable('einj','f8',('id',))  
    pinj = nc.createVariable('pinj','f8',('id','time',))        

    btilt_15  = nc.createVariable('btilt_15' ,'f8',('one',))  
    stilt_15L = nc.createVariable('stilt_15L','f8',('one',))  
    stilt_15R = nc.createVariable('stilt_15R','f8',('one',))  
    
    #data = empty(8,'O')
    #for k in range(nbi["num_srcs"]): data[k] = nb_list[k]
    #id[:]   = data 

    id[:,:] = nb_list
    time[:] = nbi["times"]

    einj[:]   = nbi["einj"]           
    pinj[:,:] = nbi["pinj"]           
    btilt_15[:] = [nbi["tilt_15"][0]]
    stilt_15L[:] = [nbi["tilt_15"][1]]
    stilt_15R[:] = [nbi["tilt_15"][2]]
    
    nc.close()

    #------------------------------------------------------------------
    # return

    return nbi

if __name__ == "__main__":

    shot= 147634 #153648 #153071 #147634
    nbi = get_nbi(shot)




