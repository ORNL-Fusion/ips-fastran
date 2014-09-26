#!/usr/bin/env Python

"""
 ======================================================================
 Last modified : Apr 2014
 post process for fastran
 JM
"""

import os,sys,glob,shutil,pickle

from numpy import *
from scipy import interpolate
import netCDF4
import Namelist
import zfdat

#######################################################################
# common data
#

define_jobs = {
    "current": [],
    }

define_options = {
    "use_cache" : ["bool" , False],
    }
                     
#######################################################################
# utils
#

def print_err():

    print '\nusage: xfastran_post.py job [--option=...]\n'
    print '       job  :'
    print '              current'
    print ''
    sys.exit(-1)

def call_external(exe,dir):

    pp = zproc.zproc()
    pp.setWorkDir(os.path.realpath(dir))
    pp.goWorkDir()
    os.system(exe)
    pp.backDir()

def integrate_s(rho,volp,ipol,r0,cur):

    nrho = len(rho)
    I = zeros(nrho)

    cur_c = node2cell(cur)
    for i in range(nrho-1):
        drho = rho[i+1]-rho[i]
        volp_c = 0.5*(volp[i+1]+volp[i])
        ipol_c = 0.5*(ipol[i+1]+ipol[i])
        I[i+1] = I[i] + cur_c[i]*(volp_c/ipol_c**2)*drho
    return I/(2.0*pi*r0)

def node2cell(node):

    nrho = len(node)
    cell = zeros(nrho-1)
    cell[0] = 0.5*(node[0]+node[1])
    for i in range(1,nrho-1):
        cell[i] = 2.0*node[i]-cell[i-1]
    return cell

def integrate_current():

    ncfile = 'fastran.nc'
    inmetric_file = 'inmetric'
    instate_file = glob.glob('*_ps.instate')[0]

    fastran = netCDF4.Dataset(ncfile,'r',format='NETCDF4') 
    ip   = fastran.variables['ip' ][:]
    ibs  = fastran.variables['ibs'][:]
    inb  = fastran.variables['inb'][:]
    irf  = fastran.variables['irf'][:]
    j_tot  = fastran.variables['j_tot'][:,:]
    j_nb  = fastran.variables['j_nb'][:,:]
    j_bs  = fastran.variables['j_bs'][:,:]
    j_rf  = fastran.variables['j_rf'][:,:]

    print "plasma current [MA]:",ip[-1]
    print "bootstrap current [MA,fraction]:",ibs[-1],ibs[-1]/ip[-1]
    print "nb current [MA,fraction]:",inb[-1],inb[-1]/ip[-1]
    print "rf current [MA,fraction]:",irf[-1],irf[-1]/ip[-1]

    inmetric = zfdat.read_formatted(inmetric_file)
    rho  = inmetric["RHO"]
    rhob = inmetric["RHOB"][0]*rho
    volp = inmetric["VOLP"]
    ipol = inmetric["IPOL"]
    nrho = len(rho)


    instate = Namelist.Namelist(instate_file)
    r0 = instate["instate"]["rmajor"][0]
   #r0   = inmetric["RMAJOR"][0]
    j_ec = instate["instate"]["j_ec"]
    j_ic = instate["instate"]["j_ic"]

    j_tot_int = integrate_s(rhob,volp,ipol,r0,j_tot[-1])
    j_bs_int  = integrate_s(rhob,volp,ipol,r0,j_bs [-1])
    j_nb_int  = integrate_s(rhob,volp,ipol,r0,j_nb [-1])
    j_rf_int  = integrate_s(rhob,volp,ipol,r0,j_rf [-1])
    j_ec_int  = integrate_s(rhob,volp,ipol,r0,j_ec)
    j_ic_int  = integrate_s(rhob,volp,ipol,r0,j_ic)

    j_tot_int = j_tot_int*(ip[-1]/j_tot_int[-1])
    j_bs_int  = j_bs_int*(ibs[-1]/j_bs_int[-1])
    j_nb_int  = j_nb_int*(inb[-1]/j_nb_int[-1])
    j_rf_int  = j_rf_int*(irf[-1]/j_rf_int[-1])
    tmp = irf[-1]/(j_ec_int[-1]+j_ic_int[-1])
    j_ec_int  = j_ec_int*tmp
    j_ic_int  = j_ic_int*tmp

    #print j_tot_int[-1]
    #print j_bs_int[-1]
    #print j_nb_int[-1]
    #print j_rf_int[-1]
    #print j_ec_int[-1]+j_ic_int[-1]
    #print j_ec_int[-1]
    #print j_ic_int[-1]

    f=open("current_integral.dat","w")
    for var in ["rho","total","bs","nb","ec","rf"]:
        f.write("%6s "%var)
    f.write("\n")
    for i in range(nrho):
        f.write("%6.3f "%rho[i])
        f.write("%6.3f "%j_tot_int[i])
        f.write("%6.3f "%j_bs_int[i])
        f.write("%6.3f "%j_nb_int[i])
        f.write("%6.3f "%j_ec_int[i])
        f.write("%6.3f "%j_ic_int[i])
        f.write("\n")
    f.close()

########################################################################
# driver
#

if __name__=="__main__":

    from optparse import OptionParser

    #-------------------------------------------------------------------
    # check input

    if len(sys.argv) < 2: print_err()

    job = sys.argv[1]

    if job not in define_jobs: print_err()

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
    # tasks

    if job=="current":

        integrate_current()


