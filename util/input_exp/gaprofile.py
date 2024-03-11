#!/usr/bin/env python3

"""
-----------------------------------------------------------------------
gaprofiles io
-----------------------------------------------------------------------
"""

import sys
import os
import pickle
import glob
from numpy import *
import scipy.optimize
from optparse import OptionParser
import idlsave
import netCDF4
from fastran.util.zinterp import zinterp

#----------------------------------------------------------------------
# UTILITY

def fname_d3data(rdir, ishot, itime, prefix, ext=''):
    return os.path.join(rdir, '%s%06d.%05d%s'%(prefix, ishot, itime, ext)) 

def get_shot_time(fname):
    ishot = int(os.path.basename(fname).split('.')[0][-6:])
    itime = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    return ishot, itime

def __screen_out(s):
    print (72*'#')
    print (s+":")
    print (72*'#')

#----------------------------------------------------------------------
#   READ GAPROFILES

__vmap = {'ne':'dne', 'te':'dte', 'ti':'dti', 'omega':'dtrot', 'nz':'dimp', 'rad':'prad'}
def read(v, shot, time, rdir='.'): 
    prefix = __vmap[v]

    if prefix == 'dne': 
       f = fname_d3data(rdir, shot, time, prefix)
       d = idlsave.read(f, verbose=False)
       x = d.ne_str.rho_dens[0]
       x2 = d.ne_str.psi_dens[0]
       y = d.ne_str.dens[0]
       yp = d.ne_str.densp[0]
       e = d.ne_str.dens_err[0]
       ind = d.ne_str.valid[0] > 0
       dx = d.ne_str.ne_data[0].rho[0][ind]
       dx2 = d.ne_str.ne_data[0].psi_norm[0][ind]
       dy = d.ne_str.ne_data[0].thom[0][ind]
       de = d.ne_str.ne_data[0].thom_err[0][ind]
       lasr = d.ne_str.ne_data[0].in_time_block[0]
       norm = d.ne_str.tnorm[0].norms[0]
       dy = norm[lasr][ind]*dy
       de = norm[lasr][ind]*de #<=========== corrected 2020.10.08

       ttime = d.ne_str.tnorm[0].times[0][lasr][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.ne_str.valid[0] < 0
       rx = d.ne_str.ne_data[0].rho[0][ind]
       ry = d.ne_str.ne_data[0].thom[0][ind]
       re = d.ne_str.ne_data[0].thom_err[0][ind]
       rhomax = d.rhomax
       yp = yp*rhomax
 
    elif prefix == 'dte':
       f = fname_d3data(rdir, shot, time, prefix)
       d = idlsave.read(f, verbose=False)
       x = d.te_str.rho_te[0]
       x2 = d.te_str.psi_te[0]
       y = d.te_str.te[0]
       yp = d.te_str.tep[0]
       e = d.te_str.te_err[0]
       ind = d.te_str.valid[0] > 0
       dx = d.te_str.te_data[0].rho[0][ind]
       dx2 = d.te_str.te_data[0].psi_norm[0][ind]
       dy = d.te_str.te_data[0].temp[0][ind]
       de = d.te_str.te_data[0].temp_err[0][ind]

       ttime = d.te_str.te_data[0].time1[0][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.te_str.valid[0] < 0
       rx = d.te_str.te_data[0].rho[0][ind]
       ry = d.te_str.te_data[0].temp[0][ind]
       re = d.te_str.te_data[0].temp_err[0][ind]
       rhomax = d.rhomax
       yp = yp*rhomax

    elif prefix == 'dti':
       f = fname_d3data(rdir, shot, time, prefix)
       d = idlsave.read(f, verbose=False)
       x = d.ti_str.rho_ti[0]
       x2 = d.ti_str.psi_ti[0]
       y = d.ti_str.ti[0]
       yp = d.ti_str.tip[0]
       e = d.ti_str.ti_err[0]
       ind = d.ti_str.valid[0] > 0
       dx = d.ti_str.ti_data[0].rho[0][ind]
       dx2 = d.ti_str.ti_data[0].psi[0][ind]
       dy = d.ti_str.ti_data[0].temp[0][ind]
       de = d.ti_str.ti_data[0].temp_err[0][ind]

       ttime = d.ti_str.ti_data[0].time[0][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.ti_str.valid[0] < 0
       rx = d.ti_str.ti_data[0].rho[0][ind]
       ry = d.ti_str.ti_data[0].temp[0][ind]
       re = d.ti_str.ti_data[0].temp_err[0][ind]
       rhomax = d.rhomax
       yp = yp*rhomax

    elif prefix == 'dtrot':
       f = fname_d3data(rdir, shot, time, prefix)
       d = idlsave.read(f, verbose=False)
       x = d.tor_rot_str.rho_tor_rot[0]
       x2 = d.tor_rot_str.psi_tor_rot[0]
       y = d.tor_rot_str.tor_rot_local[0]
       yp = [] #<================================ need to implement
       e = d.tor_rot_str.tor_rot_err[0]
       ind = d.tor_rot_str.valid[0] > 0
       dx = d.tor_rot_str.tor_rot_data[0].rho[0][ind]
       dx2 = d.tor_rot_str.tor_rot_data[0].psi_norm[0][ind]
       try:
           dy = d.tor_rot_str.tor_rot_data[0].omega_cor[0][ind]
       except AttributeError:
           print('no rotation correction')
           dy = d.tor_rot_str.tor_rot_data[0].omega[0][ind]
       de = d.tor_rot_str.tor_rot_data[0].omega_err[0][ind]

       ttime = d.tor_rot_str.tor_rot_data[0].time[0][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.tor_rot_str.valid[0] < 0
       rx = d.tor_rot_str.tor_rot_data[0].rho[0][ind]
       try:
           ry = d.tor_rot_str.tor_rot_data[0].omega_cor[0][ind]
       except AttributeError:
           print('no rotation correction')
           ry = d.tor_rot_str.tor_rot_data[0].omega[0][ind]
       re = d.tor_rot_str.tor_rot_data[0].omega_err[0][ind]
       rhomax = d.rhomax

    elif prefix == 'dimp':
       f = fname_d3data(rdir, shot, time, prefix,'_Carbon')
       d = idlsave.read(f,verbose=False)
       x = d.impdens_str.rho_imp[0]
       x2 = d.impdens_str.psi_imp[0]
       y = d.impdens_str.zdens[0]
       yp = d.impdens_str.zdensp[0]
       e = d.impdens_str.zdens_err[0]
       ind = d.impdens_str.valid[0] > -10
       dx = d.impdens_str.impdens_data[0].rho[0][ind]
       dx2 = d.impdens_str.impdens_data[0].psi[0][ind]
       dy = d.impdens_str.impdens_data[0].DATA[0][ind]
       de = d.impdens_str.impdens_data[0].data_err[0][ind]

       ttime = d.impdens_str.impdens_data[0].time[0][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.impdens_str.valid[0] < 0
       rx = d.impdens_str.impdens_data[0].rho[0][ind]
       ry = d.impdens_str.impdens_data[0].DATA[0][ind]
       re = d.impdens_str.impdens_data[0].data_err[0][ind]
       rhomax = d.rhomax
       yp = yp*rhomax

    elif prefix == 'prad':
       f = fname_d3data(rdir, shot, time, prefix)
       d = idlsave.read(f, verbose=False)
       x = d.prad_str.prad_rho[0]
       x2 = array([-1])
       y = d.prad_str.prad_prof[0]
       yp = [] #<================================ need to implement
       e = d.prad_str.prad_err[0]
       dx = array([-1])
       dx2 = array([-1])
       dy = array([-1])
       de = array([-1])
       dt = array([-1])
       rx = array([-1])
       ry = array([-1])
       re = array([-1])

    if len(rx)==0:
        rx = array([-1])
        ry = array([-1])
        re = array([-1])

    return {'x':x, 'x2':x2, 'y':y, 'yp':yp, 'e':e,
            'dx':dx, 'dx2':dx2, 'dy':dy, 'de':de, 'dt':dt, 'rx':rx, 'ry':ry, 're':re, 
            'rhomax':rhomax}

def readall(
        shot,
        vars=['ne', 'te', 'ti', 'omega', 'nz', 'rad'],
        rdir='.',
        tmin=0,
        tmax=10000,
        ifilltime=False): 

    prof = {} 
    for v in vars:
        print('reading '+__vmap[v])
        prof[v] = {} 
        files = glob.glob("%s/%s%06d.?????*"%(rdir, __vmap[v], shot))
        files.sort()
        for file in files:
           s, time = get_shot_time(file)
           if time < tmin or time > tmax: continue
           try:
               prof[v][time] = read(v, shot, time, rdir=rdir)
           except:
               print("error reading: "+file)
    prof['shot'] = shot
    for v in vars: print("number of %6s = %3d"%(v, len(prof[v])))

    if ifilltime:
        times = array([t for t in prof[vars[0]].keys()])
        print(vars[0], times)
        for v in vars[1:]:
            times = intersect1d(times, array( [t for t in prof[v].keys()] ))
            print('*', v, times)
        times.sort()
        print("-----------")
        print(times)

        for v in vars:
            for time in prof[v].keys():
                if time not in times:
                   print("remove %s, t = %d"%(v,time))

    if 'nz' in vars:
        prof['zeff'] = {}

        if not ifilltime:
            times = array( [t for t in prof['nz'].keys()] )
            times.sort()

        for time in times:
            ne = prof['ne'][time]['y'][0:101]
            nz = prof['nz'][time]['y'][0:101]
            i = where (nz > ne/6.0)[0]
            if len(i) > 0: print(i)
            for k in range(len(ne)):
                if nz[k]> ne[k]/6.0: print ('negative ion at ', k)

            zeff = 30.0*array(prof['nz'][time]['y'][0:101]) \
                       /array(prof['ne'][time]['y'][0:101])+1.0
            prof['zeff'][time] = \
                 { 'x': prof['nz'][time]['x'][0:101],
                   'y': zeff,
                   'e': zeros(101),
                   'dx': array([-1]),
                   'dy': array([-1]),
                   'de': array([-1]),
                   'dt': array([-1]),
                   'rx': [],
                   'ry': [],
                   're': []}

    return prof

#----------------------------------------------------------------------
#   MANIMUPLATE DATA

def mtanh(c, x, y, param=None):
    z = 2.*(c[0]-x)/c[1]
    pz = 1. + c[4]*z
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  ) 
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    return yfit-y

def fit_func(c, x, y, func, param=None):
    rval = func(c, x) - y
    return rval

def ts_shift(prof):
    times = [t for t in prof["te"].keys()]
    times.sort()
    print(times)

    fit_func = mtanh

    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    for time in times:

        te_x = prof["te"][time]["x"]
        te_y = prof["te"][time]["y"]
        wped_input = 0.1
        yped_input = te_x[90]
        fit_coeff_input = [1.0-0.5*wped_input, wped_input, yped_input, 0.0, te_y[0]-yped_input, 1.5, 1.5]
        fit_coeff, success = scipy.optimize.leastsq(fit_func, fit_coeff_input, args=(te_x, te_y))
        xmid = fit_coeff[0]
        xwid = fit_coeff[1]
        xsep = xmid+0.5*xwid
        print(time, xsep)

        te_s = zinterp(te_x, te_y)
        te_shift = zeros(len(te_y))

        ne_x = prof["ne"][time]["x"]
        ne_y = prof["ne"][time]["y"]
        ne_s = zinterp(ne_x, ne_y)
        ne_shift = zeros(len(ne_y))

        for k in range(len(te_y)):
            te_shift[k] = te_s(te_x[k]*xsep)
            ne_shift[k] = ne_s(ne_x[k]*xsep)

        prof["te"][time]["y0"] = prof["te"][time]["y"] 
        prof["ne"][time]["y0"] = prof["ne"][time]["y"] 

        prof["te"][time]["y"] = te_shift 
        prof["ne"][time]["y"] = ne_shift

def time_avg(
    prof,
    time0,
    dtavg,
    vars = ['ne', 'te', 'ti', 'omega', 'nz'],
    include_rad=False,
    include_zeff=False):

    if include_rad: vars = vars+['rad']
    if include_zeff: vars = vars+['zeff']
    nrho = 101
    rho = arange(nrho)/float(nrho-1)
    prof_avg = {} 
    prof_avg['shot'] = prof['shot']
    for v in vars:
        times = array(prof[v].keys())
        times.sort()
        times = times[logical_and(times>time0-dtavg,times<time0+dtavg)]
        ntimes = len(times)
        if ntimes < 1: 
            print("no profile [%5s] in t = [%d,%d]"%(v,time0-dtavg,time0+dtavg))
            print("...exit")
            sys.exit(-1)
        p = zeros((ntimes,nrho))
        for k in range(ntimes):
            p[k][:nrho] = prof[v][times[k]]['y'][:nrho]
        prof_avg[v] = {}
        prof_avg[v][time0] = {}
        prof_avg[v][time0]['x'] = rho
        prof_avg[v][time0]['y'] = average(p,axis=0)
        prof_avg[v][time0]['e'] = std(p,axis=0)
        prof_avg[v][time0]['dx'] = array([])
        prof_avg[v][time0]['dy'] = array([])
        prof_avg[v][time0]['de'] = array([])
        for k in range(ntimes):
            for d in ['dx','dy','de']:
                prof_avg[v][time0][d] = \
                    append(prof_avg[v][time0][d],prof[v][times[k]][d])
        print("profile average: [%5s]"%v,time0,dtavg,times)
    return prof_avg

#----------------------------------------------------------------------
#   NETCDF

def wrt_netcdf(shot, prof, outdir='.', ncfile=''):
    if ncfile:
        print('append to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'r+', format='NETCDF4_CLASSIC')
    else:
        ncfile = 'p%06d.nc'%ishot
        print('dump to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')

    vars = [] 
    for var in prof.keys():
        if var in ['ne', 'te', 'ti', 'omega', 'nz', 'rad']: 
            vars.append(var)

    ndata = {}
    for var in vars:
        ndata[var] = {}
        for key in ["dx", "dy", "de", "dt", "rx", "ry", "re"]:
            ndata[var][key] = 0 
            for time in prof[var]:
                ndata[var][key] = max( len(prof[var][time][key]), ndata[var][key] )
    for key in ndata:
        print(key, ndata[key])

    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    nc_rho = {}
    for var in vars:
        nc.createDimension('rho_%s'%var , nrho)
        nc_rho[var]  = nc.createVariable('rho_%s'%var ,'f8', ('rho_%s'%var, ))
        nc_rho[var][:] = rho[:]

    nc_time = {}
    for var in vars:
        times = sorted(prof[var].keys())
        nc.createDimension('time_%s'%var, len(times))
        nc_time[var] = nc.createVariable('time_%s'%var, 'i4', ('time_%s'%var, ))
        nc_time[var][:] = times[:]

    #for var in vars:
    #    nc.createDimension('ndata_%s'%var, len(ndata[var]["dx"]))
        
    nc_vars = {}
    for var in vars:
        nc_vars[var] = nc.createVariable(var, 'f8', ('time_%s'%var, 'rho_%s'%var))      

    #for var in vars:
    #    for key in ["dx", "dy", "de", "dt"]:
    #        nc_vars[var+"_"+key] = nc.createVariable(var+"_"+key,'f8',('time_%s'%var,"ndata_"+var))      

    for var in vars:
        for k, time in enumerate(nc_time[var][:]):
            try:
                nc_vars[var][k,:] = prof[var][time]['y'][:nrho]
            except:
                print('no data', var)
                nc_vars[var][k,:] = zeros(nrho)
            #for key in ["dx", "dy", "de", "dt"]:
            #    n = len(nc_vars[var+"_"+key][k,:])
            #    out = array(prof[var][time][key][:])
            #    out.resize(n)
            #    nc_vars[var+"_"+key][k,:] = out[:] #prof[var][time][key][:]

    nc.close()

    from Namelist import Namelist

#----------------------------------------------------------------------
# TEST

define_jobs = {
    "dump"   : [],
    "plot"   : [],
#    "fit"    : [],
    "check"  : [],
    }

define_options = {
    "rad" : ["bool" , False],
    "ts_shift" : ["bool" , False],
    }
                     
def print_err():
    print('\nusage: ped-profile.py job shot tmin tmax prof_dir [--option=...]\n')
    print('       job  :')
    print('              dump, plot, check')
    print('')
    raise Exception("input argument error")


def input_args():
    parser = OptionParser()
    for option in define_options:
        type = define_options[option][0]
        default = define_options[option][1]
        if type == "bool": 
            parser.add_option("--"+option,
                action="store_true",dest=option,default=False)
        else:
            parser.add_option("--"+option,
                action="store",type=type,dest=option,default=default)
    (options,args) = parser.parse_args(sys.argv[1:])
    return options

if __name__=="__main__":
    if len(sys.argv) < 2: print_err()
    job = sys.argv[1]
    if job not in define_jobs: print_err()

    shot = int(sys.argv[2])
    tmin = int(sys.argv[3])
    tmax = int(sys.argv[4])
    prof_dir = sys.argv[5]

    options = input_args()

    if job=="dump":
        if options.rad:
            prof = readall(shot, vars=['ne', 'te', 'ti', 'omega', 'nz', 'rad'], rdir=prof_dir)
        else:
            prof = readall(shot, vars=['ne', 'te', 'ti', 'omega', 'nz'], rdir=prof_dir)
        if options.ts_shift:
            ts_shift(prof)
        pickle.dump(prof,open("prof%06d_%s"%(shot,os.path.basename(prof_dir)), "wb"))
        wrt_netcdf(shot, prof, outdir='.')
