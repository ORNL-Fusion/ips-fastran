import os
import sys
import shutil
import glob
import multiprocessing
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

from collections import deque

def run_ips_wrapper(arg):
    print ('RUN', arg)

    try:
        run_ips(arg)
        print ('RETURN NORMAL')
        return 0
    except:
        print ('RETURN EXCEPT')
        return -1 

def run_ips(arg):
    k, f_sim_config, header, data, cleanup  = arg[0], arg[1], arg[2], arg[3], arg[4]

    cfgFile = os.path.realpath("run%05d.config"%k)
    rundir = os.path.realpath("run%05d"%k)
    logFile = "ipslog.%05d"%k
    platformFile = os.path.realpath("local.conf")

    if not os.path.exists(rundir): os.makedirs(rundir)
    if not os.path.exists("output"): os.makedirs("output")
    output_dir = os.path.realpath("output")

    sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

    for i, key in enumerate(header.split()):
        print (i,key)
        comp, vname, vtype  = key.split(':') 
        d = data.split() 
        if vtype == 'str': val = str(d[i])
        else: val = eval(vtype+"("+d[i]+")")
        if comp == '': sim[vname] = val
        else: sim[comp][vname] = val 

    sim["SIM_ROOT"] = rundir

    sim.write(open(cfgFile, "wb"))

    t0 = timer.time()
    olddir = os.getcwd()
    os.chdir(rundir)

    try:
        do_create_runspace  =  True
        do_run_setup        =  True
        do_run              =  True
        compset_list        =  []
        debug               =  None
        ftb                 =  None
        verbose_debug       =  None
        cmd_nodes           =  0
        cmd_ppn             =  0

        print(('START %5d'%k))
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

        del fwk
    except:
        print(('FAILED',rundir))

    os.chdir(olddir)
    t1 = timer.time()
    print(('FINISHED %5d %6.3f'%(k, t1-t0)))
    sys.stdout.flush()

    dir_state = sim['PLASMA_STATE_WORK_DIR']
    for filename in glob.glob(os.path.join(dir_state, "*.*")):
        shutil.copy(filename, output_dir)
    if cleanup:
        shutil.rmtree(rundir)

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

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(("*Starting IPS, massive serial", rank, size, os.uname()[1]))
    sys.stdout.flush()
    
    if len(sys.argv) < 3: 
        print ("usage: ips_massive_serial.py inscan sim_config ppn")
        sys.exit()

    f_inscan = sys.argv[1]
    f_sim_config = sys.argv[2]
    ppn = int(sys.argv[3])
    cleanup = int(sys.argv[4])
      
    with open(f_inscan,"r") as f: data = f.readlines()
    nsim = len(data)-1

    wrt_localconf("local.conf")

    scan_list = deque()
    for k in range(rank, nsim, size): scan_list.append( [k, f_sim_config, data[0], data[k+1], cleanup, rank] )
    nsim = len(scan_list)
    print(('NSIM = %d at RANK = %d'%(nsim, rank)))

    procs = []

    for k in range(0, ppn):
        if len(scan_list)>0:
            scan = scan_list.popleft()
            p = multiprocessing.Process(target= run_ips_wrapper, args=(scan,)) 
            procs.append(p)
            p.start() 

    while (len(scan_list)>0):

         timer.sleep(5)
         print(('PROCS LENTH=', len(procs)))
         sys.stdout.flush()

         procs_next = []

         for k in range(len(procs)):
             p = procs[k]
             if p.is_alive():
                procs_next.append(p)
             else:
                p.join()
                if len(scan_list)>0:
                    scan = scan_list.popleft()
                    p1 = multiprocessing.Process(target=run_ips_wrapper, args=(scan,)) 
                    procs_next.append(p1)
                    p1.start()

         procs = procs_next
         
    comm.Barrier()

    print ('END ALL')

if __name__ == "__main__":
    main()

