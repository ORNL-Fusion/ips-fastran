#!/usr/bin/python

import os
import sys
import shutil

import multiprocessing
import multiprocessing.pool

from mpi4py import MPI

import time as timer

import ips
from configobj import ConfigObj

import numpy
import scipy
import netCDF4
import plasmastate
import Namelist
import subprocess

class NoDaemonProcess(multiprocessing.Process):
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def run_ips(arg):

    k, f_sim_config, header, data, cleanup  = arg[0], arg[1], arg[2], arg[3], arg[4]

    cfgFile = os.path.realpath("run%05d.config"%k)
    rundir = os.path.realpath("run%05d"%k)
    logFile = "ipslog.%05d"%k
    platformFile = os.path.realpath("local.conf")

    if not os.path.exists(rundir): os.makedirs(rundir)

    sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

    for i, key in enumerate(header.split()):
        comp, vname, vtype  = key.split(':') 
        d = data.split() 
        if vtype == 'str': val = str(d[i])
        else: val = eval(vtype+"("+d[i]+")")
        if comp == '': sim[vname] = val
        else: sim[comp][vname] = val 

    sim["SIM_ROOT"] = rundir

    sim.write(open(cfgFile,"w"))

    t0 = timer.time()
    olddir = os.getcwd()
    os.chdir(rundir)

    do_create_runspace  =  True
    do_run_setup        =  True
    do_run              =  True
    compset_list        =  []
    debug               =  None
    ftb                 =  None
    verbose_debug       =  None
    cmd_nodes           =  0
    cmd_ppn             =  0

    print 'START %5d'%k
    print platformFile
    sys.stdout.flush()

    fwk = ips.Framework(
        do_create_runspace,
        do_run_setup,
        do_run,
        [cfgFile],
        logFile,
        platformFile,
        compset_list,
        debug,
        ftb,
        verbose_debug,
        cmd_nodes,
        cmd_ppn)

    fwk.run()

    os.chdir(olddir)
    t1 = timer.time()
    print 'FINISHED %5d %6.3f'%(k,t1-t0)
    sys.stdout.flush()

    t0 = timer.time()
    dir_state = sim['PLASMA_STATE_WORK_DIR']
    shutil.copytree(dir_state, 'out%05d'%k)
    if cleanup:
        shutil.rmtree(rundir)
    t1 = timer.time()
    print 'DIR_STATE %s %6.3f'%(dir_state,t1-t0)
    sys.stdout.flush()

    return 0

def wrt_localconf(fname="local.conf"):

    s = \
"""HOST = edison
MPIRUN = eval
NODE_DETECTION = manual
TOTAL_PROCS = 1
NODES = 1
PROCS_PER_NODE = 1
CORES_PER_NODE = 1
SOCKETS_PER_NODE = 1
NODE_ALLOCATION_MODE = SHARED
USE_ACCURATE_NODES = ON
"""
    f=open(fname,"w")
    f.write(s)


if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print "*Starting IPS, massive serial",rank
    sys.stdout.flush()
    
    if len(sys.argv) < 3: 
        print "usage: ips_massive_serial.py inscan sim_config ppn"
        sys.exit()

    f_inscan = sys.argv[1]
    f_sim_config = sys.argv[2]
    ppn = int(sys.argv[3])

    try:
        cleanup = int(sys.argv[4])
    except:
        cleanup = 0
      
    with open(f_inscan,"r") as f: data = f.readlines()
    nsim = len(data)-1

    wrt_localconf("local.conf")

    processes_list = []
    for k in range(rank,nsim,size): processes_list.append([k,f_sim_config,data[0],data[k+1],cleanup] )
    print processes_list
    print 'NSIM = %d at RANK = %d'%(len(processes_list),rank)
    pool = NoDaemonProcessPool(processes = ppn)
    result = pool.map_async(run_ips, processes_list)
    output = result.get()

    comm.Barrier()

    print 'END ALL'
