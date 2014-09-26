#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 collect efit geqdsks 
 last modified Jun, 2013
 ----------------------------------------------------------------------
"""

import os,sys,pickle,glob
from numpy import *
from scipy import interpolate
import Namelist
import pmds

import netCDF4
import zefitutil,zefitdata
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

def get_shot_time(fname):
    shot = int(os.path.basename(fname).split('.')[0][-6:])
    time = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    return (shot,time)

def get_efit(shot,tmin=2000,tmax=6000,outdir='.',efitdir='.'):

    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    nbdry_max = 201

    gfiles = glob.glob(os.path.join(efitdir,'g%06d.?????'%shot))
    gfiles.sort()

    times = []
    geqs = []

    for gfile in gfiles:

        time = get_shot_time(gfile)[-1]
        if time<tmin or time>tmax: continue

        print 'processing: %s'%gfile
        times.append(time)
       #geqs.append(zefitutil.readg(gfile))
        geqs.append(zefitdata.efitdata(gfile,"ps"))

    ntimes = len(times)
    ip = zeros(ntimes)
    b0 = zeros(ntimes)
    r0 = zeros(ntimes)
    rmajor = zeros(ntimes)
    aminor = zeros(ntimes)
    kappa = zeros(ntimes)
    delta = zeros(ntimes)
    r = zeros((ntimes,nrho))
    a = zeros((ntimes,nrho))
    jtot = zeros((ntimes,nrho))
    p_eq = zeros((ntimes,nrho))
    q_eq = zeros((ntimes,nrho))
    nbdry = zeros(ntimes,'i4')
    nlim  = zeros(ntimes,'i4')
    rbdry = zeros((ntimes,nbdry_max))
    zbdry = zeros((ntimes,nbdry_max))
    rlim = zeros((ntimes,nbdry_max))
    zlim = zeros((ntimes,nbdry_max))
    
    for k in range(ntimes):
        geq = geqs[k]
        ip[k] = geq["cpasma"]
        b0[k] = geq["bcentr"]
        r0[k] = geq["rzero"]
        rmajor[k] = geq["ps"]["Rmajor_mean"][-1]
        aminor[k] = geq["ps"]["rMinor_mean"][-1]
        kappa[k] = geq["ps"]["elong"][-1]
        delta[k] = geq["ps"]["triang"][-1]
        jtot[k,:] = geq["ps"]["jdotb"]/b0[k]
        p_eq[k,:] = geq["ps"]["P_eq"]
        q_eq[k,:] = geq["ps"]["q_eq"]
        r[k,:]    = geq["ps"]["Rmajor_mean"]
        a[k,:]    = geq["ps"]["rMinor_mean"]

        bdry = zefitutil.getbdry(geq,nskip=1)
        nbdry[k] = bdry["nbdry"]
        nlim[k] = bdry["nlim"]
        rbdry[k,0:nbdry[k]] = bdry["rbdry"] #[0:nbdry[k]]
        zbdry[k,0:nbdry[k]] = bdry["zbdry"]
        rlim[k,0:nlim[k]] = bdry["rlim"]
        zlim[k,0:nlim[k]] = bdry["zlim"]

    #------------------------------------------------------------------
    # dump to netcdf

    print 'dump to '+outdir+'/geq_%06d.nc'%shot

    nc = netCDF4.Dataset(outdir+'/geq_%06d.nc'%shot,'w'
             ,format=xfastran_env.netcdf_form)
    
    dim_rho  = nc.createDimension('rho' , nrho)
    dim_time = nc.createDimension('time', ntimes)
    dim_bdry = nc.createDimension('bdry', nbdry_max)
    
    var_rho    = nc.createVariable('rho','f8',('rho',))
    var_time   = nc.createVariable('time','f8',('time',))

    var_ip     = nc.createVariable('ip','f8',('time',))
    var_b0     = nc.createVariable('b0','f8',('time',))
    var_r0     = nc.createVariable('r0','f8',('time',))
    var_rmajor = nc.createVariable('rmajor','f8',('time',))
    var_aminor = nc.createVariable('aminor','f8',('time',))
    var_kappa  = nc.createVariable('kappa','f8',('time',))
    var_delta  = nc.createVariable('delta','f8',('time',))
    var_jtot   = nc.createVariable('jtot','f8',('time','rho'))
    var_p_eq   = nc.createVariable('p_eq','f8',('time','rho'))
    var_q_eq   = nc.createVariable('q_eq','f8',('time','rho'))
    var_r      = nc.createVariable('r','f8',('time','rho'))
    var_a      = nc.createVariable('a','f8',('time','rho'))

    var_nbdry  = nc.createVariable('nbdry','i4',('time',))
    var_nlim   = nc.createVariable('nlim','i4',('time',))
    var_rbdry  = nc.createVariable('rbdry','f8',('time','bdry'))
    var_zbdry  = nc.createVariable('zbdry','f8',('time','bdry'))
    var_rlim   = nc.createVariable('rlim','f8',('time','bdry'))
    var_zlim   = nc.createVariable('zlim','f8',('time','bdry'))

    var_rho[:]    = rho
    var_time[:]   = times
    var_ip[:]     = ip
    var_b0[:]     = b0
    var_r0[:]     = r0
    var_rmajor[:] = rmajor
    var_aminor[:] = aminor
    var_kappa[:]  = kappa
    var_delta[:]  = delta
    var_jtot[:,:] = jtot
    var_p_eq[:,:] = p_eq
    var_q_eq[:,:] = q_eq
    var_r[:,:] = r
    var_a[:,:] = a

    var_nbdry[:]  = nbdry
    var_nlim[:]   = nlim
    var_rbdry[:,:] = rbdry
    var_zbdry[:,:] = zbdry
    var_rlim[:,:] = rlim
    var_zlim[:,:] = zlim

    nc.close()

    #------------------------------------------------------------------
    # return

    return 0

if __name__ == "__main__":

    shot= 153646
    efitdir = '/u/parkjm/013/20132204/153646/kin02'
   #tmp = get_efit(shot,3700,4700,outdir='.',efitdir=efitdir)
    tmp = get_efit(shot,3700,3800,outdir='.',efitdir=efitdir)




