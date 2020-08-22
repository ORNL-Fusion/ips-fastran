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
import ipsutil

class ips_massive_serial(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0):
        self.clean_after = int(getattr(self, "CLEAN_AFTER", "0"))
        self.time_out = int(getattr(self, "TIME_OUT", "3600000"))

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

        tmp_xfs = getattr(self, "TMPXFS", "")
        # check tmp xfs filesystem is setup correctly
        if tmp_xfs:
            if os.system(f'findmnt -nt xfs -T {tmp_xfs}') != 0:
                raise Exception(f"TMPXFS is set but this is either not running in a shifter container "
                                f"or {tmp_xfs} is not mounted as a temporary xfs file")

        f_sim_config = getattr(self, "SIMULATION", "")
        sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

        f_inscan = getattr(self, "INSCAN", "inscan")
        with open(f_inscan,"r") as f: inscan = f.readlines()
        nsim = len(inscan)-1
        header = inscan[0]

        dask_nodes = int(getattr(self, "DASK_NODES", "1"))
        task_ppn = int(getattr(self, "TASK_PPN", ""))

        rank = int(getattr(self, "RANK", "0"))
        size = int(getattr(self, "SIZE", "1"))
        print(('NSIM = %d at RANK = %d'%(nsim, rank)))

        try: 
            pwd = services.get_config_param("PWD")
        except:
            pwd = os.environ["PWD"]

        #--- add to pool
        pool = services.create_task_pool('pool')
        cwd = services.get_working_dir()

        if tmp_xfs:
            tmp_xfs_rank = os.path.join(tmp_xfs, str(rank)) # make unique folder per rank to work in the tmp xfs directory
            # copy files needed by ips_massive_serial.py
            ipsutil.copyFiles(".", getattr(self, "INPUT_FILES"), tmp_xfs_rank)
            # copy files needed for fastran_modeleq
            ipsutil.copyFiles(os.path.join(pwd, "input"), "*", os.path.join(tmp_xfs_rank, "input"))
            ipsutil.copyFiles(".", "local.conf", tmp_xfs_rank)
            tmp_xfs_dir_summary = os.path.realpath(os.path.join(tmp_xfs_rank, "SUMMARY"))
            if not os.path.exists(tmp_xfs_dir_summary): os.makedirs(tmp_xfs_dir_summary)

        tasks = {}
        for k in range(rank, nsim, size):
            rundir = os.path.realpath(os.path.join(tmp_xfs_rank, "run%05d"%k)) if tmp_xfs else os.path.realpath("run%05d"%k)
            logfile = os.path.realpath(os.path.join(tmp_xfs_rank, "ipslog.%05d"%k)) if tmp_xfs else "ipslog.%05d"%k
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

            sim["PWD"] = tmp_xfs_rank if tmp_xfs else pwd
            sim["RUN_ID"] = "run%05d"%k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(tmp_xfs_rank if tmp_xfs else cwd, "run%05d.out"%k)
            sim["USE_PORTAL"] = "False"
            driver = sim['PORTS']['DRIVER']['IMPLEMENTATION'] 
            sim[driver]["SUMMARY"] = tmp_xfs_dir_summary if tmp_xfs else dir_summary
            sim.write(open(os.path.join(tmp_xfs_rank, "run%05d.config"%k), "wb") if tmp_xfs else open("run%05d.config"%k, "wb"))

            ips_bin = os.path.join(self.BIN_PATH, self.BIN) 
            ips_bin += " --config=run%05d.config"%k
            ips_bin += " --log=ips_%05d.log"%k
            ips_bin += " --platform=local.conf"
            if self.clean_after: 
                ips_bin += "\nrm -rf "+rundir
                ips_bin += "\nrm -f "+rundir+".zip"

            with open(os.path.join(rundir, "ips_bin.sh"), "w") as f: f.write(ips_bin)

            services.add_task(
                'pool', 
                'task_'+str(k), 
                dask_nodes, 
                tmp_xfs_rank if tmp_xfs else cwd,
                "sh", 
                os.path.join(rundir, "ips_bin.sh"),
                logfile=logfile,
                timeout=self.time_out,
                task_ppn=task_ppn)

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

