#!/usr/bin/env Python

"""
 ======================================================================
 Last modified : Aug 2013
 driver for fastran/ips
 JM
"""

import os,sys,glob,shutil,pickle,re

from numpy import *
from scipy import interpolate
import Namelist
import netCDF4
import zechdata,znbidata,zgaprofile,zgeqdsk,zefitdata,zefitutil
import zd3dutil,znubeaminput

#######################################################################
# common data
#

verid = "*** xfastran.py ver1.0/Aug2013\n"
dir_fastran  = os.environ["FASTRAN_ROOT"]
dir_template = os.path.join(dir_fastran,'template')

#######################################################################
# submit jobs on venus and lohan
#
def submitjob_venus(inputdir):

    inputdir = os.path.realpath(inputdir)
    
    config = "fastran_scenario.config"

    f=open(dir_template+"/submitjob_venus","r")
    lines = f.readlines()
    f.close()

    command =  "Python $IPS_ROOT/bin/ips --config=%s --platform=venus.conf"%config

    f=open("submitjob_venus","w")
    for k,line in enumerate(lines):
        if line[0]=="#":
            f.write(line)
    f.write(command)
    f.close()

    os.system("qsub submitjob_venus")

#######################################################################
# driver
#

if __name__=="__main__":

    pass
