"""
 -----------------------------------------------------------------------
 massive serial run
 -----------------------------------------------------------------------
"""

import os
import shutil
import numpy as np
import time
from configobj import ConfigObj
from component import Component

class ips_massive_serial(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0):
        self.timeout = float(getattr(self, "TIME_OUT", "3600"))

    def step(self, timeStamp=0):
        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- setup
        wrt_localconf("local.conf")

        dir_summary = getattr(self, "SUMMARY", "SUMMARY")
        dir_summary = os.path.realpath(dir_summary)
        if not os.path.exists(dir_summary): os.makedirs(dir_summary)
         
        f_sim_config = getattr(self, "SIMULATION", "")
        sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

        f_inscan = getattr(self, "INSCAN", "inscan")
        with open(f_inscan,"r") as f: inscan = f.readlines()
        nsim = len(inscan)-1
        header = inscan[0]

        dask_nodes = int(getattr(self, "DASK_NODES", "-1"))
        if dask_nodes < 0:
           raise Exception("DASK_NODES undefined")

        rank = int(getattr(self, "RANK", "0"))
        size = int(getattr(self, "SIZE", "1"))
        print(('NSIM = %d at RANK = %d'%(nsim, rank)))

        clean_after = int(getattr(self, "CLEAN_AFTER", "0"))

        #--- add to pool
        pool = services.create_task_pool('pool')
        cwd = services.get_working_dir()
        tasks = {}
        for k in range(rank, nsim, size): 
            rundir = os.path.realpath("run%05d"%k)
            logfile = "ipslog.%05d"%k
            if not os.path.exists(rundir): os.makedirs(rundir)

            data = inscan[k+1]
            for i, key in enumerate(header.split()):
                print (i, key)
                comp, vname, vtype = key.split(':') 
                d = data.split() 
                if vtype == 'str': val = str(d[i])
                else: val = eval(vtype+"("+d[i]+")")
                if comp == '': sim[vname] = val
                else: sim[comp][vname] = val 

            sim["PWD"] = os.environ["PWD"]
            sim["RUN_ID"] = "run%05d"%k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(cwd, "run%05d.out"%k)
            sim["USE_PORTAL"] = "False"
            driver = sim['PORTS']['DRIVER']['IMPLEMENTATION'] 
            sim[driver]["SUMMARY"] = dir_summary
            sim.write(open("run%05d.config"%k, "wb"))

            ips_bin = os.path.join(self.BIN_PATH, self.BIN) 
            ips_bin += " --config=run%05d.config"%k
            ips_bin += " --log=ips_%05d.log"%k
            ips_bin += " --platform=local.conf"
            if clean_after: 
                ips_bin += "\nrm -rf "+rundir
                ips_bin += "\nrm -f "+rundir+".zip"

            with open(os.path.join(rundir, "ips_bin.sh"), "w") as f: f.write(ips_bin)

            services.add_task(
                'pool', 
                'task_'+str(k), 
                1, 
                cwd,
                "sh", 
                os.path.join(rundir, "ips_bin.sh"),
                logfile=logfile,
                timeout=-1)

        #--- run
        ret_val = services.submit_tasks('pool', use_dask=True, dask_nodes=dask_nodes)
        print('ret_val = ', ret_val)
        exit_status = services.get_finished_tasks('pool')
        print(exit_status)

        #--- update plasma state files
        services.update_plasma_state()

        #--- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0):
        pass  

def wrt_localconf(fname="local.conf"):
    s = \
"""HOST = local
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
    f=open(fname, "w")
    f.write(s)
    f.close()

