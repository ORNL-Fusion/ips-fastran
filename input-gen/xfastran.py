#!/usr/bin/env python

"""
 ======================================================================
 input generator for fastran/ips
"""

import os
import sys
import shutil
import cPickle
from numpy import *

import Namelist
import netCDF4
import zechdata
import znbidata
import zgaprofile
import zgeqdsk
import xfastran_input
import xfastran_env

#----------------------------------------------------------------------
# common data
#

define_jobs = {
    "input"  : [],
    "run"    : [],
    "plot"   : [],
    "check"  : []
    }

define_options = {
    "use_cache" : ["bool" , False],
    "shot" : ["int" , -1],
    "tmin" : ["int" , -1],
    "tmax" : ["int" , -1],
    "prof_dir" : ["string" , ''],
    "efit_dir" : ["string" , ''],
    "efit_id" : ["string" , ''],
    "time_plasma_shape" : ["int" , -1],
    "instate" : ["bool", False],
    "innubeam": ["bool", False],
    "intoray" : ["bool", False],
    "ts_shift" : ["bool", False],
    }
                     
#----------------------------------------------------------------------
# utils
#

def print_err():

    print '\nusage: xfastran.py job [--option=...]\n'
    print '       job  :'
    print '              input, run, plot, check'
    print ''
    raise Exception("input argument error")

def call_external(exe,dir):

    pp = zproc.zproc()
    pp.setWorkDir(os.path.realpath(dir))
    pp.goWorkDir()
    os.system(exe)
    pp.backDir()

#----------------------------------------------------------------------
# read parameters
#

def generate_input(input,infiles):

    # read inscenario

    shot = input["in_d3d"]["shot"][0]
    tmin = input["in_d3d"]["tmin"][0]
    tmax = input["in_d3d"]["tmax"][0]

    prof_dir = input["in_d3d"]["prof_dir"][0]
    efit_dir = input["in_d3d"]["efit_dir"][0]
    input_dir = input["in_d3d"]["input_dir"][0]
    data_dir = input["in_d3d"]["data_dir"][0]

    time_plasma_shape = input["in_d3d"]["time_plasma_shape"][0]

    if not os.path.exists(data_dir): os.makedirs(data_dir)
    if not os.path.exists(input_dir): os.makedirs(input_dir)

    #-------------------------------------------------------------------
    # collect data

    if not options.use_cache:
       try:
           zechdata.get_ech(shot,tmin,tmax,outdir=data_dir)
           iech = True
       except:
           iech = False
       znbidata.get_nbi(shot,tmin,tmax,outdir=data_dir)
       if prof_dir:
           prf=zgaprofile.readall(shot,rdir=prof_dir,tmin=tmin,tmax=tmax)
           if options.ts_shift:
               zgaprofile.ts_shift(prf)
           zgaprofile.wrt_netcdf(shot,prf,outdir=data_dir)
           cPickle.dump(prf,open("prof_%06d"%shot,"wb"))
       if efit_dir:
           zgeqdsk.get_efit(shot,tmin,tmax,outdir=data_dir,efitdir=efit_dir,efit_id=options.efit_id)

    #-------------------------------------------------------------------
    # write intoray

    if "intoray" in infiles and iech:
        print 'writing intoray'
        xfastran_input.wrt_intoray(data_dir+'/ech_%06d.nc'%shot
            ,shot,tmin,tmax,input_dir)
        shutil.copyfile("intoray","intoray_%06d"%shot)

    #-------------------------------------------------------------------
    # write innubeam

    if "innubeam" in infiles:
        print 'writing innubeam'
        xfastran_input.wrt_innubeam(data_dir+'/nbi_%06d.nc'%shot
            ,shot,tmin,tmax,input_dir
            ,Db=0.0)
        shutil.copyfile("innubeam","innubeam_%06d"%shot)
        print 'end'

    #-------------------------------------------------------------------
    # write instate

    if "instate" in infiles:
        print 'writing instate'
        xfastran_input.wrt_instate(
            data_dir+'/prf_%06d.nc'%shot,
            data_dir+'/geq_%06d.nc'%shot,
            shot,tmin,tmax,time_plasma_shape,input_dir)
        shutil.copyfile("instate","instate_%06d"%shot)

    #-------------------------------------------------------------------
    # write infastran

    if "infastran" in infiles:
        print 'writing infastran'
        xfastran_input.wrt_infastran(input,input_dir)

    #-------------------------------------------------------------------
    # write ingeqdsk

    if "ingeqdsk" in infiles:
        print 'writing ingeqdsk'
        xfastran_input.wrt_ingeqdsk(input_dir)

    #-------------------------------------------------------------------
    # write configuration file

    #print 'writing configuration file'
    #xfastran_input.wrt_config(input,'.')

    #-------------------------------------------------------------------
    # clean-up temporary files

    for file in ["bdry.dat","geqdsk_state","inps","log_tmp"]:
        if os.path.exists(file): 
            os.remove(file) 

def generate_input_tdep(input,infiles):

    #-------------------------------------------------------------------
    # read inscenario

    shot = input["in_d3d"]["shot"][0]
    tmin = input["in_d3d"]["tmin"][0]
    tmax = input["in_d3d"]["tmax"][0]

    prof_dir = input["in_d3d"]["prof_dir"][0]
    efit_dir = input["in_d3d"]["efit_dir"][0]
    input_dir = input["in_d3d"]["input_dir"][0]
    data_dir = input["in_d3d"]["data_dir"][0]

    time_plasma_shape = input["in_d3d"]["time_plasma_shape"][0]

    if not os.path.exists(data_dir): os.makedirs(data_dir)
    if not os.path.exists(input_dir): os.makedirs(input_dir)

    #-------------------------------------------------------------------
    # collect data

    if not options.use_cache:
       try:
           zechdata.get_ech(shot,tmin,tmax,outdir=data_dir)
           iech = True
       except:
           iech = False
       znbidata.get_nbi(shot,tmin,tmax,outdir=data_dir)
       if prof_dir:
           prf=zgaprofile.readall(shot,rdir=prof_dir,tmin=tmin,tmax=tmax)
           zgaprofile.wrt_netcdf(shot,prf,outdir=data_dir)
       if efit_dir:
           zgeqdsk.get_efit(shot,tmin,tmax,outdir=data_dir,efitdir=efit_dir)

    #-------------------------------------------------------------------
    # write intoray

    if "intoray" in infiles and iech:
        print 'writing intoray'
        xfastran_input.wrt_intoray(data_dir+'/ech_%06d.nc'%shot
            ,shot,tmin,tmax,input_dir)

    #-------------------------------------------------------------------
    # write innubeam

    if "innubeam" in infiles:
        print 'writing innubeam'
        xfastran_input.wrt_innubeam(data_dir+'/nbi_%06d.nc'%shot
            ,shot,tmin,tmax,input_dir
            ,Db=0.0)
        print 'end'

    #-------------------------------------------------------------------
    # write instate

    ncfile_prf = data_dir+'/prf_%06d.nc'%shot
    nc_prf = netCDF4.Dataset(ncfile_prf,'r',format='NETCDF4')
    time_prf = nc_prf.variables["time"][:]

    print time_prf

    if "instate" in infiles:
        for k in range(1,len(time_prf)-1):
            tmin = time_prf[k-1]+1.0
            tmax = time_prf[k+1]-1.0 
            print 'writing instate', time_prf[k]
            xfastran_input.wrt_instate(
                data_dir+'/prf_%06d.nc'%shot,
                data_dir+'/geq_%06d.nc'%shot,
                shot,tmin,tmax,time_plasma_shape,input_dir,fn_instate="i%06d.%05d"%(shot,int(time_prf[k])))

    #-------------------------------------------------------------------
    # write infastran

    if "infastran" in infiles:
        print 'writing infastran'
        xfastran_input.wrt_infastran(input,input_dir)

    #-------------------------------------------------------------------
    # write ingeqdsk

    if "ingeqdsk" in infiles:
        print 'writing ingeqdsk'
        xfastran_input.wrt_ingeqdsk(input_dir)

    #-------------------------------------------------------------------
    # write configuration file

    #print 'writing configuration file'
    #xfastran_input.wrt_config(input,'.')

    #-------------------------------------------------------------------
    # clean-up temporary files

    for file in ["bdry.dat","geqdsk_state","inps","log_tmp"]:
        if os.path.exists(file): 
            os.remove(file) 


########################################################################
# driver
#

if __name__=="__main__":

    from optparse import OptionParser

    #-------------------------------------------------------------------
    # check input

    if len(sys.argv) < 2: print_err()

    job = sys.argv[1]

    if job not in ["input","run","check"]: print_err()

    #-------------------------------------------------------------------
    # parse options

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

    #-------------------------------------------------------------------
    # read inscenario

    # input = Namelist.Namelist("inscenario")

    #-------------------------------------------------------------------
    # run

    if job=="input":

       input = Namelist.Namelist()

       input["in_d3d"]["shot"] = [options.shot]
       input["in_d3d"]["tmin"] = [options.tmin]
       input["in_d3d"]["tmax"] = [options.tmax]

       input["in_d3d"]["prof_dir" ] = [options.prof_dir]
       input["in_d3d"]["efit_dir" ] = [options.efit_dir]
       input["in_d3d"]["input_dir"] = ['.']
       input["in_d3d"]["data_dir" ] = ['.']

       input["in_d3d"]["time_plasma_shape"] =  [options.time_plasma_shape]


      #infiles = ["innubeam","intoray","instate"]
       infiles = []
       if options.innubeam: infiles += ["innubeam"]
       if options.intoray: infiles += ["intoray"]
       if options.instate: infiles += ["instate"]
       print infiles

       generate_input(input,infiles)
     # generate_input_tdep(input,infiles)

    elif job=="run":

       print 'please use job batch system'
