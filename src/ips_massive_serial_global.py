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

class ips_massive_serial_global(Component):

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
        dir_summary = getattr(self, "SUMMARY", "SUMMARY")
        dir_summary = os.path.realpath(dir_summary)
        if not os.path.exists(dir_summary): os.makedirs(dir_summary)
         
        f_sim_config = getattr(self, "SIMULATION", "")
        sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

        f_machine_config = getattr(self, "MACHINE", "")

        ntasks = int(getattr(self, "NTASKS", "1"))
        dask_nodes = 1 
        task_ppn = 1

        print('NTASKS = %d'%ntasks)

        #--- add to pool
        pool = services.create_task_pool('pool')
        cwd = services.get_working_dir()
        tasks = {}
        for k in range(ntasks): 
            rundir = os.path.realpath("run%05d"%k)
            logfile = "ipslog.%05d"%k
            if not os.path.exists(rundir): os.makedirs(rundir)

            sim["PWD"] = os.environ["PWD"]
            sim["RUN_ID"] = "run%05d"%k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(cwd, "run%05d.out"%k)
            sim["USE_PORTAL"] = "False"

            driver = sim['PORTS']['DRIVER']['IMPLEMENTATION'] 
            sim[driver]["SUMMARY"] = dir_summary#+"%05d"%k
            sim[driver]["RANK"] = k
            sim[driver]["SIZE"] = ntasks 
            sim.write(open("run%05d.config"%k, "wb"))

            ips_bin = os.path.join(self.BIN_PATH, self.BIN) 
            ips_bin += " --config=run%05d.config"%k
            ips_bin += " --log=ips_%05d.log"%k
            ips_bin += " --platform=%s"%f_machine_config
            if self.clean_after: 
                ips_bin += "\nrm -rf "+rundir
                ips_bin += "\nrm -f "+rundir+f_machine_config

            with open(os.path.join(rundir, "ips_bin.sh"), "w") as f: f.write(ips_bin)

            services.add_task(
                'pool', 
                'task_'+str(k), 
                dask_nodes, 
                cwd,
                "sh", 
                os.path.join(rundir, "ips_bin.sh"),
                logfile=logfile,
                timeout=self.time_out,
                task_ppn=task_ppn
            )

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
