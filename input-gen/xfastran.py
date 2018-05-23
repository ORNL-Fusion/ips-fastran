#!/usr/bin/env Python

"""
 ======================================================================
 Last modified : Aug 2013
 driver for fastran/ips
 JM
"""

import os,sys,glob,shutil,pickle

from numpy import *
from scipy import interpolate
import Namelist
import netCDF4
import zechdata,znbidata,zgaprofile,zgeqdsk,zefitdata,zefitutil
import zd3dutil,znubeaminput
import xfastran_input,xfastran_submitjob
import xfastran_env

#######################################################################
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
    }
                     
#######################################################################
# utils
#

def print_err():

    print '\nusage: xfastran.py job [--option=...]\n'
    print '       job  :'
    print '              input, run, plot, check'
    print ''
    sys.exit(-1)

def call_external(exe,dir):

    pp = zproc.zproc()
    pp.setWorkDir(os.path.realpath(dir))
    pp.goWorkDir()
    os.system(exe)
    pp.backDir()

########################################################################
# read parameters
#

def generate_input(input):

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
       prf=zgaprofile.readall(shot,rdir=prof_dir)
       zgaprofile.wrt_netcdf(shot,prf,outdir=data_dir)
       zgeqdsk.get_efit(shot,tmin,tmax,outdir=data_dir,efitdir=efit_dir)

    #-------------------------------------------------------------------
    # write intoray

    if iech:
        print 'writing intoray'
        xfastran_input.wrt_intoray(data_dir+'/ech_%06d.nc'%shot
            ,shot,tmin,tmax,input_dir)

    #-------------------------------------------------------------------
    # write innubeam

    print 'writing innubeam'
    xfastran_input.wrt_innubeam(data_dir+'/nbi_%06d.nc'%shot
        ,shot,tmin,tmax,input_dir
        ,Db=0.0)

    #-------------------------------------------------------------------
    # write instate

    print 'writing instate'
    xfastran_input.wrt_instate(
        data_dir+'/prf_%06d.nc'%shot,
        data_dir+'/geq_%06d.nc'%shot,
        shot,tmin,tmax,time_plasma_shape,input_dir)

    #-------------------------------------------------------------------
    # write infastran

    print 'writing infastran'
    xfastran_input.wrt_infastran(input,input_dir)

    #-------------------------------------------------------------------
    # write ingeqdsk

    print 'writing ingeqdsk'
    xfastran_input.wrt_ingeqdsk(input_dir)

    #-------------------------------------------------------------------
    # write configuration file

    print 'writing configuration file'
    xfastran_input.wrt_config(input,'.')

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

    input = Namelist.Namelist("inscenario")

    #-------------------------------------------------------------------
    # run

    if job=="input":

       generate_input(input)

    elif job=="run":

       xfastran_submitjob.submitjob_venus('.')
