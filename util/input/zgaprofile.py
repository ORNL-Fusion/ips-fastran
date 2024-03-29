#!/usr/bin/env Python

"""
-----------------------------------------------------------------------
Manipulate gaprofiles
last modified : May 2010
JM

Requirements:
* numpy
* idlsave
* Namelist

Supports:
* generate inone (tdem, time slice, time average, ...)
* generate ufile for transp
* new fit in time (box car average, smooth, linear fit, spline fit, ....)
* elm average
* sawtooth average 
-----------------------------------------------------------------------
"""

import sys
import os
import cPickle
import glob
from collections import OrderedDict
from numpy import *
import scipy.optimize
from optparse import OptionParser
import idlsave
import netCDF4
from zinterp import zinterp
import xfastran_env

#----------------------------------------------------------------------
# UTILITY
#

def fname_d3data(rdir,shot,time,prefix,ext=''):

    return rdir+'/'+'%s%06d.%05d'%(prefix,shot,time)+ext 

def get_shot_time(fname):

    shot = int(os.path.basename(fname).split('.')[0][-6:])
    #time = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    tmp = os.path.basename(fname).split('.')[-1].split('_')[0]
    k = 0
    for v in tmp:
        if v=='0': k=k+1
        else: break
    time = int(tmp[k:])
    return shot,time

def __screen_out(s):

    print 72*'#'
    print s+":"
    print 72*'#'

#----------------------------------------------------------------------
#   READ GAPROFILES
#

__vmap = {'ne':'dne','te':'dte','ti':'dti','omega':'dtrot','nz':'dimp'
    ,'rad':'prad'}

def read(v,shot,time,rdir='.'): 

    prefix = __vmap[v]

    if prefix == 'dne': 

       f = fname_d3data(rdir,shot,time,prefix)
       d = idlsave.read(f,verbose=False)
       x = d.ne_str.rho_dens[0]
       x2 = d.ne_str.psi_dens[0]
       y = d.ne_str.dens[0]
       e = d.ne_str.dens_err[0]
       ind = d.ne_str.valid[0] > 0
       dx = d.ne_str.ne_data[0].rho[0][ind]
       dx2 = d.ne_str.ne_data[0].psi_norm[0][ind]
       dy = d.ne_str.ne_data[0].thom[0][ind]
       de = d.ne_str.ne_data[0].thom_err[0][ind]
       lasr = d.ne_str.ne_data[0].in_time_block[0]
       norm = d.ne_str.tnorm[0].norms[0]
       dy = norm[lasr][ind]*dy

       ttime = d.ne_str.tnorm[0].times[0][lasr][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.ne_str.valid[0] < 0
       rx = d.ne_str.ne_data[0].rho[0][ind]
       ry = d.ne_str.ne_data[0].thom[0][ind]
       re = d.ne_str.ne_data[0].thom_err[0][ind]
      
    elif prefix == 'dte':

       f = fname_d3data(rdir,shot,time,prefix)
       d = idlsave.read(f,verbose=False)
       x = d.te_str.rho_te[0]
       x2 = d.te_str.psi_te[0]
       y = d.te_str.te[0]
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

    elif prefix == 'dti':

       f = fname_d3data(rdir,shot,time,prefix)
       d = idlsave.read(f,verbose=False)
       x = d.ti_str.rho_ti[0]
       x2 = d.ti_str.psi_ti[0]
       y = d.ti_str.ti[0]
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

    elif prefix == 'dtrot':

       f = fname_d3data(rdir,shot,time,prefix)
       d = idlsave.read(f,verbose=False)
       x = d.tor_rot_str.rho_tor_rot[0]
       x2 = d.tor_rot_str.psi_tor_rot[0]
       y = d.tor_rot_str.tor_rot_local[0]
       e = d.tor_rot_str.tor_rot_err[0]
       ind = d.tor_rot_str.valid[0] > 0
       dx = d.tor_rot_str.tor_rot_data[0].rho[0][ind]
       dx2 = d.tor_rot_str.tor_rot_data[0].psi_norm[0][ind]
       try:
           dy = d.tor_rot_str.tor_rot_data[0].omega_cor[0][ind]
       except AttributeError:
           print 'no rotation correction'
           dy = d.tor_rot_str.tor_rot_data[0].omega[0][ind]
       de = d.tor_rot_str.tor_rot_data[0].omega_err[0][ind]

       ttime = d.tor_rot_str.tor_rot_data[0].time[0][ind]
       dt = array([ int(ttime[k]) for k in range(len(ttime)) ])

       ind = d.tor_rot_str.valid[0] < 0
       rx = d.tor_rot_str.tor_rot_data[0].rho[0][ind]
       try:
           ry = d.tor_rot_str.tor_rot_data[0].omega_cor[0][ind]
       except AttributeError:
           print 'no rotation correction'
           ry = d.tor_rot_str.tor_rot_data[0].omega[0][ind]
       re = d.tor_rot_str.tor_rot_data[0].omega_err[0][ind]

    elif prefix == 'dimp':

       f = fname_d3data(rdir,shot,time,prefix,'_Carbon')
       d = idlsave.read(f,verbose=False)
       x = d.impdens_str.rho_imp[0]
       x2 = d.impdens_str.psi_imp[0]
       y = d.impdens_str.zdens[0]
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

    elif prefix == 'prad':

       f = fname_d3data(rdir,shot,time,prefix)
       d = idlsave.read(f,verbose=False)
       x = d.prad_str.prad_rho[0]
       x2 = array([-1])
       y = d.prad_str.prad_prof[0]
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

    return {'x':x,'x2':x2,'y':y,'e':e,'dx':dx,'dx2':dx2,'dy':dy,'de':de,'dt':dt,
            'rx':rx,'ry':ry,'re':re}

def readall(shot,vars=['ne','te','ti','omega','nz','rad'],rdir='.',
        tmin=0,tmax=10000): 
#def readall(shot,vars=['ne','te','ti','omega','nz'],rdir='.',
#        tmin=0,tmax=10000): 

    print rdir
    prof = OrderedDict()
    for v in vars:
        print 'reading '+__vmap[v]
        prof[v] = OrderedDict()
        files = glob.glob("%s/%s%06d.?????*"%(rdir,__vmap[v],shot))
        files.sort()
        for file in files:
           s, time = get_shot_time(file)
           if time < tmin or time > tmax: continue
           prof[v][time] = read(v,shot,time,rdir=rdir)
    prof['shot'] = shot
    for v in vars: print "number of %6s = %3d"%(v,len(prof[v]))

    if 'nz' in vars:
        prof['zeff'] = {}
        times = prof['nz'].keys()
        times.sort()
        for time in times:
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
#

def mtanh(c,x,y,param=None):
    z = 2.*(c[0]-x)/c[1]
    pz = 1. + c[4]*z
    mth = 0.5*( ( pz + 1.0 )*tanh(z) + pz - 1.0  ) 
    yfit = 0.5*( (c[2]-c[3])*mth + c[2]+c[3] )
    return yfit-y

def fit_func(c,x,y,func,param=None):
    rval = func(c,x)-y
    return rval

def ts_shift(prof):

    times = prof["te"].keys()
    times.sort()
    print times

    fit_func = mtanh

    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    for time in times:

        te_x = prof["te"][time]["x"]
        te_y = prof["te"][time]["y"]
        wped_input = 0.1
        yped_input = te_x[90]
        fit_coeff_input = [1.0-0.5*wped_input,wped_input,yped_input,0.0, te_y[0]-yped_input, 1.5, 1.5]
        fit_coeff, success = scipy.optimize.leastsq(fit_func,fit_coeff_input,args=(te_x,te_y))
        xmid = fit_coeff[0]
        xwid = fit_coeff[1]
        xsep = xmid+0.5*xwid
        print time, xsep

        te_s = zinterp(te_x,te_y)
        te_shift = zeros(len(te_y))

        ne_x = prof["ne"][time]["x"]
        ne_y = prof["ne"][time]["y"]
        ne_s = zinterp(ne_x,ne_y)
        ne_shift = zeros(len(ne_y))

        for k in range(len(te_y)):
            te_shift[k] = te_s(te_x[k]*xsep)
            ne_shift[k] = ne_s(ne_x[k]*xsep)

        prof["te"][time]["y0"] = prof["te"][time]["y"] 
        prof["ne"][time]["y0"] = prof["ne"][time]["y"] 

        prof["te"][time]["y"] = te_shift 
        prof["ne"][time]["y"] = ne_shift

def time_avg(prof,time0,dtavg,
    vars = ['ne','te','ti','omega','nz'],
    include_rad=False,include_zeff=False):

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
            print "no profile [%5s] in t = [%d,%d]"%(v,time0-dtavg,time0+dtavg)
            print "...exit"
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
        print "profile average: [%5s]"%v,time0,dtavg,times
    return prof_avg

def linear(t,t0,t1,y0,y1):

    nj = len(y0)
    y = zeros(nj)

    t = float(t)
    t0 = float(t0)
    t1 = float(t1)

    for j in range(nj):
        f0 = y0[j]
        f1 = y1[j]
        y[j] = f0*(t-t1)/(t0-t1) + f1*(t-t0)/(t1-t0)

    return y

def time_interpolation(p,include_rad=False): #,tmin=None,tmax=None):

    time = {}
    vars = ['ne','te','ti','omega','nz']
    if include_rad: vars = vars+['rad']
    for v in vars:
        tmp = p[v].keys()
        tmp.sort()
        time[v] = array(tmp)
       #print v,time[v]
    
    if len(time['ne']) != len(time['te']) :
       print 'ne, te time slices different...exit',len(time['ne']),len(time['te'])
    #   sys.exit()

    if include_rad:
       if len(time['ne']) != len(time['rad']) :
           print 'ne, rad time slices different...exit',len(time['ne']),len(time['rad'])
           sys.exit()
    
    timeall = array(time['ne'])
    #if tmin == None: tmin = timeall[0]
    #if tmax == None: tmax = timeall[-1]
    #id = timeall>=tmin and timeall<=tmax
    #timeall = timeall[id]

    for v in ['te','ti','omega','nz']:
        for t in time['ne']:
            if t not in time[v]:
                t0 = time[v][where( time[v]<t)][-1]
                t1 = time[v][where( time[v]>t)][0]
                print '%s: mssing t = %04d...interpolation %04d %04d'%(v,t,t0,t1)
                y0 = p[v][t0]['y']
                y1 = p[v][t1]['y']
                x = p[v][t0]['x']
                y = linear(t,t0,t1,y0,y1)
                p[v][t] = {'x':x,'y':y}
    return p

#----------------------------------------------------------------------
#   PLOT UTILS
#

import matplotlib
matplotlib.rcParams['axes.linewidth']=0.5
#matplotlib.rcParams['pdf.fonttype']=42
from matplotlib.ticker import MultipleLocator
import pylab as p

from matplotlib.backends.backend_pdf import PdfPages

def plot(shot,rdir,prof=None,yrange=None
    ,icommon_time=None,include_rad=False):

    if prof==None:
        __screen_out("READ GAPROFILES")
        prof = readall(shot,rdir=rdir)

    vars = ['ne','te','ti','omega','nz']
    if include_rad: vars = vars + ['rad']
    times = array(prof[vars[0]].keys())
    if icommon_time: 
        for v in vars[1:]:
            times = intersect1d_nu(times,array(prof[v].keys()))
    times.sort()
    print 'time: \n', times

    if yrange==None:
        yrange = {}
        for var in vars:
            ymax = 0
            for time in times:
                if prof[var].has_key(time):
                    ymax = max([ymax,max(prof[var][time]['y'])])
            yrange[var] = array([0.0,ymax])
        yrange["omega"] = 1.0e-5*yrange["omega"]
        for var in vars:
            if yrange[var][1] < 1.0:
                yrange[var][1]=ceil(12.0*yrange[var][1])/10.0
            else:
                yrange[var][1]=ceil(1.2*yrange[var][1])
            print var,yrange[var]

    pdf = PdfPages('profile_%06d.pdf'%shot)

    __screen_out("PLOT")
    for time in times:
        print 't = ',time
        plot_page([[prof,[time]]],
            mode=['data','dataerr','removed','removederr'],
            colors=['k'],
            filename="tmp%04d.eps"%time,yrange=yrange,
            include_rad=include_rad)

        pdf.savefig()

    pdf.close()

def plot_page(
    cases,
    mode = [],
    colors=['b','r','g','m','y'],
    filename="profile_page.pdf",
    yrange={"ne":[0,8.],"te":[0,8.],"ti":[0,8.],"omega":[0,1.5],"nz":[0,0.3],"rad":[0,0.3]},
    include_rad=None):

    p.figure(figsize=(8.5,8.5))
  # p.suptitle('t = %d'%time, fontsize=10)
  # p.figtext(0.8,0.05,'userid:parkjm',size='small')
    offset = 0.15
    kplot = 0
    for case in cases:
        prof = case[0]
        for time in case[1]:
            if kplot < len(colors) : color = colors[kplot]
            else: color = 'k'
            p.figtext(0.12+offset*kplot,0.95,'%d@%d'%(prof['shot'],time),color=color)
            aplot(321,prof,time,'ne',mode,color,yrange["ne"][1],1.0,'ne')
            aplot(323,prof,time,'te',mode,color,yrange["te"][1],1.0,'Te')
            aplot(325,prof,time,'ti',mode,color,yrange["ti"][1],1.0,'Ti')
            aplot(324,prof,time,'omega',mode,color,yrange["omega"][1],1.0e-5,'Omega')
            aplot(322,prof,time,'nz',mode,color,yrange["nz"][1],1.0,'nz')
            if include_rad:
                aplot(326,prof,time,'rad',['fiterr'],color,yrange["rad"][1],1.0,'rad')
            aplot(326,prof,time,'te',mode,color,yrange["te"][1],1.0,'Te')
            aplot(326,prof,time,'ti',mode,color,yrange["ti"][1],1.0,'Ti')
            kplot = kplot+1

def aplot(pos,prof,time,v,mode,color,y_max,y_scale,title):

    if not time in prof[v]: return

    ax = p.subplot(pos)

    n_tick_major = 5
    n_tick_minor = 4
    markersize = 2

    x_majL = MultipleLocator(0.20)
    x_minL = MultipleLocator(0.04)
    y_majL = MultipleLocator(y_max/n_tick_major)
    y_minL = MultipleLocator(y_max/n_tick_major/4)

    x  = prof[v][time]['x' ]
    y  = prof[v][time]['y' ]*y_scale
    e  = prof[v][time]['e' ]*y_scale
    dx = prof[v][time]['dx']
    dy = prof[v][time]['dy']*y_scale
    de = prof[v][time]['de']*y_scale
    if 'removed' in mode:
        rx = prof[v][time]['rx']
        ry = prof[v][time]['ry']*y_scale
        re = prof[v][time]['re']*y_scale

    ax.plot(x,y,color+'-',lw=1)

    if 'fiterr' in mode:
        ax.errorbar(x[::4],y[::4],e[::4],fmt=color+'o',ms=0)

    if 'data' in mode: 
        if dx[0] != -1:
            ax.plot(dx,dy,color+'o',
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=markersize)

    if 'dataerr' in mode:
        if dx[0] != -1:
            ax.errorbar(dx,dy,de,fmt=color+'o',
                markerfacecolor=color,
                markeredgecolor=color,
                ms=markersize)

    if 'removed' in mode:
        if rx[0] != -1:
            ax.plot(rx,ry,'w'+'o',
                markerfacecolor='w',
                markeredgecolor=color,
                ms=markersize)

    if 'removederr' in mode:
        if rx[0] != -1:
            ax.errorbar(rx,ry,re,fmt=color+'o',
                markerfacecolor='w',
                markeredgecolor=color,
                ms=markersize)

    ax.xaxis.set_major_locator(x_majL)
    ax.xaxis.set_minor_locator(x_minL)
    ax.yaxis.set_major_locator(y_majL)
    ax.yaxis.set_minor_locator(y_minL)
    ax.axis([0,1,0,y_max])
    p.xticks(fontsize=10)
    p.yticks(fontsize=10)
    ax.set_title(title,fontsize=10)

def tplot(profs,v,mode,colors,y_max,y_scale):

    filename="time_page.eps"
    p.figure(figsize=(8.5,8.5))
    offset = 0.15
    kplot = 0

    for prof in profs:

        color = colors[kplot]

        times  = prof[v].keys()
        times.sort()
    
        p.figtext(0.12+offset*kplot,0.95,'%d'%(prof['shot']),color=color)
    
        for i in [0,20,40,60,80,100]:
    
            title = "%3d"%i
            pos = 320+i/20+1
            ax = p.subplot(pos)
            y = []
            e = []
            for time in times:
                y.append (prof[v][time]['y'][i]*y_scale)
                e.append (prof[v][time]['e'][i]*y_scale)
        
            ax.plot(times,y,color+'-',lw=1)
    
            ax.axis([0,times[-1],0,y_max])
            p.xticks(fontsize=10)
            p.yticks(fontsize=10)
            ax.set_title(title,fontsize=10)

        kplot = kplot+1

    p.savefig(filename)
    p.clf()

def plot_example_average(shot,rdir,time0,dtavg):

    __screen_out("READ GAPROFILES")
    shot = 136345
    prof = readall(shot,rdir=rdir)

    __screen_out("PLOT")

    time0 = 3000
    dtavg = 500
    prof_0 = time_avg(prof,time0,dtavg)
    time1 = 3000
    dtavg = 500
    prof_1 = time_avg(prof,time1,dtavg)
    plot_page( [ [prof_0,[time0]], [prof_1,[time1]] ],
               mode=['fiterr','data'],
               colors=['b','r']
              )

#----------------------------------------------------------------------
#   NETCDF
#

def wrt_netcdf(shot,prof,outdir='.'):

    print 'dump to '+outdir+'/prf_%06d.nc'%shot

    nc = netCDF4.Dataset(outdir+'/prf_%06d.nc'%shot,'w'
             ,format=xfastran_env.netcdf_form)
    
    times = prof['ne'].keys()
    times.sort()
    ntime = len(times)
    nrho = 101
    rho = arange(nrho)/(nrho-1.0)

    dim_time = nc.createDimension('time', ntime)
    dim_rho  = nc.createDimension('rho' , nrho)
    
    nc_time = nc.createVariable('time','f8',('time',))
    nc_rho  = nc.createVariable('rho' ,'f8',('rho',))

    vars=['ne','te','ti','omega','nz','rad']
    nc_vars = {}
    for var in vars:
        nc_vars[var] = nc.createVariable(var,'f8',('time','rho'))      
    
    nc_time[:] = times[:]
    nc_rho[:] = rho[:]

    for var in vars:
        for k,time in enumerate(times):
            try:
                nc_vars[var][k,:] = prof[var][time]['y'][:nrho]
            except:
                print 'no data', var 
                nc_vars[var][k,:] = zeros(nrho)

    nc.close()

#----------------------------------------------------------------------
#   TEST
#

define_jobs = {
    "dump"   : [],
    "plot"   : [],
    "fit"    : [],
    "check"  : [],
    }

define_options = {
    "rad" : ["bool" , True],
    "ts_shift" : ["bool" , False],
    }
                     
def print_err():

    print '\nusage: ped-profile.py job shot tmin tmax prof_dir [--option=...]\n'
    print '       job  :'
    print '              dump, plot,fit, check'
    print ''
    sys.exit(-1)


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

        prof = readall(shot,rdir=prof_dir)
        if options.ts_shift:
            ts_shift(prof)
        cPickle.dump(prof,open("prof%06d_%s"%(shot,os.path.basename(prof_dir)),"wb"))
        wrt_netcdf(shot,prof,outdir='.')

