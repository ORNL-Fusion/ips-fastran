#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 collect efit geqdsks 
 last modified Jun, 2013
 ----------------------------------------------------------------------
"""

import os,sys,pickle,glob,shutil
from numpy import *
from scipy import interpolate
import Namelist
import pmds

import netCDF4
import efit_eqdsk,plasmastate
from zinterp import zinterp

import xfastran_env

def cell2node(cell):

    nrho = len(cell)+1
    node = zeros(nrho)
    node[0] = cell[0]
    for i in range(1,nrho-1):
        node[i] = 0.5*(cell[i-1]+cell[i])
    node[-1] = cell[-1]
    return node

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

def get_efit(shot,tmin=2000,tmax=6000,outdir='.',efitdir='.',efit_id=''):

    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    nbdry_max = 201

    gfiles = glob.glob(os.path.join(efitdir,'g%06d.?????%s'%(shot,efit_id)))
    gfiles.sort()

    times = []
    for gfile in gfiles:
        time = get_shot_time(gfile)[-1]
        if time<tmin or time>tmax: continue
        times.append(time)

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
    jpar = zeros((ntimes,nrho))
    jtor = zeros((ntimes,nrho))
    p_eq = zeros((ntimes,nrho))
    q_eq = zeros((ntimes,nrho))
    nbdry = zeros(ntimes,'i4')
    nlim  = zeros(ntimes,'i4')
    rbdry = zeros((ntimes,nbdry_max))
    zbdry = zeros((ntimes,nbdry_max))
    rlim = zeros((ntimes,nbdry_max))
    zlim = zeros((ntimes,nbdry_max))
    r_out = zeros((ntimes,nrho))
    r_in = zeros((ntimes,nrho))

    
    for k in range(ntimes):

        gfile = os.path.join(efitdir,'g%06d.%05d%s'%(shot,times[k],efit_id))
        print 'processing',gfile

        bin_fluxsurf = os.path.join(os.environ["FLUXSURF_BIN_PATH"],os.environ["FLUXSURF_BIN_NAME"])
        print bin_fluxsurf
        os.system(bin_fluxsurf+" %s"%gfile)
        current = loadtxt("current.dat").transpose()
        jpar[k,:] = current[1]
        jtor[k,:] = current[2]

        g  = efit_eqdsk.readg(gfile)
        ps = plasmastate.plasmastate('ips',1)
        ps.init_from_geqdsk(gfile)
   
        ip[k] = g["cpasma"]
        b0[k] = g["bcentr"]
        r0[k] = g["rzero"]
        rmajor[k] = ps["Rmajor_mean"][-1]
        aminor[k] = ps["rMinor_mean"][-1]
        kappa[k]  = ps["elong"][-1]
        delta[k]  = ps["triang"][-1]
        jtot[k,:] = ps["jdotb"][:]/b0[k]
        p_eq[k,:] = ps["P_eq"][:]
        q_eq[k,:] = ps["q_eq"][:]
        r[k,:]    = ps["Rmajor_mean"][:]
        a[k,:]    = ps["rMinor_mean"][:]
        bdry = efit_eqdsk.getbdry(g,nskip=1)
        nbdry[k] = bdry["nbdry"]
        nlim[k] = bdry["nlim"]
        rbdry[k,0:nbdry[k]] = bdry["rbdry"] 
        zbdry[k,0:nbdry[k]] = bdry["zbdry"]
        rlim[k,0:nlim[k]] = bdry["rlim"]
        zlim[k,0:nlim[k]] = bdry["zlim"]

        r_out[k,:] = ps["R_midp_out"][:]
        r_in[k,:] = ps["R_midp_in"][:]

        psi = ps["psipol"][:]/ps["psipol"][-1] 
      # rho = sqrt(ps["phit"][:]/ps["phit"][-1])
        rhob = (ps["phit"][-1]/pi/abs(b0[k]))**0.5 
        rhopsi = zinterp(psi,rho)
        ipol = ps["g_eq"][:]/(r0[k]*b0[k])
        ipol = abs(ipol)

        #ipol_tmp = g["fpol"]/(r0[k]*b0[k])
        #psin = arange(len(ipol_tmp))/(len(ipol_tmp)-1.)
        #print len(psin)
        #ipol_spl = zinterp(psin,ipol_tmp)
        #ipol = ipol_spl(psi)
        #print ipol

        volp = 4.0*pi*pi*rho*rhob*r0[k]/ipol/(r0[k]*r0[k]*ps["gr2i"][:])
        g11 = volp*ps["grho2"][:]*rhob**2
        g22 = r0[k]*volp/(4.0*pi*pi)*ps["grho2r2i"][:]*rhob**2
        g33 = r0[k]*r0[k]*ps["gr2i"][:]

#        psipol = 2.0*pi*ps["psipol"][:]
#       #psipol = ( g["ssibry"] - g["ssimag"] ) * arange(nrho)/(nrho-1.0) +  g["ssimag"]
#        drho = rhob/(nrho-1.0)
#        jpar = zeros(nrho)
#        for i in range(1,nrho-1):
#           jpar[i] = (g22[i]+g22[i+1])*(psipol[i+1]-psipol[i])-(g22[i-1]+g22[i])*(psipol[i]-psipol[i-1])
#           jpar[i] *= 5.*ipol[i]**2/volp[i]*0.5/drho**2
#        jpar[-1] = 2.0*jpar[-2]-jpar[-3]
#        jpar[0] = 2.0*jpar[1]-jpar[2]
#        jtot[k,:] = jpar*1.0e6
       

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
    var_jpar   = nc.createVariable('jpar','f8',('time','rho'))
    var_jtor   = nc.createVariable('jtor','f8',('time','rho'))
    var_p_eq   = nc.createVariable('p_eq','f8',('time','rho'))
    var_q_eq   = nc.createVariable('q_eq','f8',('time','rho'))
    var_r      = nc.createVariable('r','f8',('time','rho'))
    var_a      = nc.createVariable('a','f8',('time','rho'))

    var_r_out   = nc.createVariable('r_out','f8',('time','rho'))
    var_r_in   = nc.createVariable('r_in','f8',('time','rho'))

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
    var_jpar[:,:] = jpar
    var_jtor[:,:] = jtor
    var_p_eq[:,:] = p_eq
    var_q_eq[:,:] = q_eq
    var_r[:,:] = r
    var_a[:,:] = a

    var_r_out[:,:] = r_out
    var_r_in[:,:] = r_in

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




