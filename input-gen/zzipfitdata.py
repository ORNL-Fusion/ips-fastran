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

def boxsmooth(y,n):
    ys = zeros(len(y))
    for i in range(n,nlen(y)-n):
        tmp = 0
        for k in range(2*n+1):
            tmp+=y[i-n-k]
        tmp = tmp/(2*n+1.0)
        ys[i] = ys
    return ys

def get_zipfit(shot,tmin=2000,tmax=6000,dt=20,outdir='.'):

    #------------------------------------------------------------------
    # connect mds and open nb tree

    print 'connecting to mds server'
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    base_tree = 'electrons'

    print 'open electrons tree'
    pmds.mdsopen(base_tree,shot)

    #------------------------------------------------------------------
    # get data

    print 'get data'

    node = '.profile_fits.zipfit'
    vname = 'edensfit'

    print node+':'+vname
    ne = pmds.mdsvalue(node+':'+vname) #ne[0] : first time slice, ....
    print ne[:,0]

    time = pmds.mdsvalue('DIM_OF(\\%s,0)'%vname)
    rho  = pmds.mdsvalue('DIM_OF(\\%s,1)'%vname)
    print time
    print rho
   
    prefix+shot.time+suffix

    dne+"%06d.%05d"%(shot,time)


    # time_mds = pmds.mdsvalue('DIM_OF(\\PINJ)')
    # pinj_mds = {}
    # einj_mds = {}
    # for src in nb_list:
    #     pinj_mds[src] = pmds.mdsvalue('\\PINJ_'+src)
    #     einj_mds[src] = pmds.mdsvalue('.NB'+src+':NBVAC_SCALAR')

    # try:
    #     btilt_15  = pmds.mdsvalue('.NB15L.OANB:BLPTCH_CAD')
    #     stilt_15L = pmds.mdsvalue('.NB15L.OANB:SRCPTCH')/60.0
    #     stilt_15R = pmds.mdsvalue('.NB15R.OANB:SRCPTCH')/60.0
    # except:
    #     btilt_15  = 0.0
    #     stilt_15L = 0.0
    #     stilt_15R = 0.0

    # print len(time_mds), len(pinj_mds['30L']), einj_mds['30L']
    # print time_mds[0],time_mds[1],time_mds[1]-time_mds[0]
    # print time_mds[-1],time_mds[-2],time_mds[-1]-time_mds[-2]
    # print btilt_15,stilt_15L,stilt_15R
    #  
    # #------------------------------------------------------------------
    # # process data

    # times = arange(tmin,tmax,dt)
    # num_times = len(times)

    # pinj = zeros((len(nb_list),num_times))
    # einj = zeros(len(nb_list))
    # nbsrc = len(nb_list)*[""]

    # for k, src in enumerate(nb_list):

    #     for i in range(num_times):
    #         ind = (time_mds >= times[i]-0.5*dt ) & (time_mds <= times[i]+0.5*dt )
    #         pinj[k,i] = average(pinj_mds[src][ind])
    #         if pinj[k,i] <= 3.5e4: pinj[k,i]=0.0

    #     einj[k] = einj_mds[src] 

    #     nbsrc[k] = src

    # #------------------------------------------------------------------
    # # disconnect mds

    # pmds.mdsdisconnect()

    # #------------------------------------------------------------------
    # # check plot

    # p = zplot.zplot("nbi_%06d.pdf"%shot,size_x=8.0,size_y=4.0)

    # for k, src in enumerate(nb_list):
    #     zplot.plot_s(111,1.0e-3*times,1.0e-6*pinj[k,:]
    #          ,[0.0,7.0,1.0],[0.0,3.0,1.0]
    #         ,iline='line',isym='',xlab='time',ylab='pinj_%s'%src)
    #     p.savefig()

    # p.close()


    # #------------------------------------------------------------------
    # # output

    # nbi = { 
    #     "num_srcs"    : len(nb_list),  
    #     "nbsrc"       : nbsrc,
    #     "num_times"   : num_times,  
    #     "times"       : times,      # [msec]
    #     "pinj"        : pinj,       # [W]
    #     "einj"        : einj,       # [eV]
    #     "tilt_15"     : [btilt_15,stilt_15L,stilt_15R] #[deg]
    # }

    # #------------------------------------------------------------------
    # # dump to netcdf

    # nc = netCDF4.Dataset(outdir+'/nbi_%06d.nc'%shot,'w',format='NETCDF4')
    # 
    # dim_id   = nc.createDimension('id'  , nbi["num_srcs"])
    # dim_time = nc.createDimension('time', nbi["num_times"])
    # dim_one  = nc.createDimension('one' , 1)
    # 
    # id   = nc.createVariable('id',str,('id',))
    # time = nc.createVariable('time','f8',('time',))

    # einj = nc.createVariable('einj','f8',('id',))  
    # pinj = nc.createVariable('pinj','f8',('id','time',))        

    # btilt_15  = nc.createVariable('btilt_15' ,'f8',('one',))  
    # stilt_15L = nc.createVariable('stilt_15L','f8',('one',))  
    # stilt_15R = nc.createVariable('stilt_15R','f8',('one',))  
    # 
    # data = empty(8,'O')
    # for k in range(nbi["num_srcs"]): data[k] = nb_list[k]
    # id[:]   = data 
    # time[:] = nbi["times"]

    # einj[:]   = nbi["einj"]           
    # pinj[:,:] = nbi["pinj"]           
    # btilt_15[:] = [nbi["tilt_15"][0]]
    # stilt_15L[:] = [nbi["tilt_15"][1]]
    # stilt_15R[:] = [nbi["tilt_15"][2]]
    # 
    # nc.close()

    #------------------------------------------------------------------
    # return

    # return nbi

if __name__ == "__main__":

    shot= 153646
    nbi = get_zipfit(shot)




