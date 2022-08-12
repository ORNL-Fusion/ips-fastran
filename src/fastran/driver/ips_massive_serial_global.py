"""
 -----------------------------------------------------------------------
 massive serial run
 -----------------------------------------------------------------------
"""
import os
import shutil
import numpy as np
import time
import subprocess
from ipsframework.configobj import ConfigObj
from ipsframework import Component

class ips_massive_serial_global(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        self.clean_after = int(getattr(self, "CLEAN_AFTER", "0"))
        self.time_out = int(getattr(self, "TIME_OUT", "3600000"))

    def step(self, timeid=0):
        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

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
        print('NTASKS = %d'%ntasks)

        #--- generate simulation config file for each node
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
            sim["USE_PORTAL"] = "True"
            sim['PARENT_PORTAL_RUNID'] = self.services.get_config_param("PORTAL_RUNID")

            driver = sim['PORTS']['DRIVER']['IMPLEMENTATION']
            sim[driver]["SUMMARY"] = dir_summary
            sim[driver]["RANK"] = k
            sim[driver]["SIZE"] = ntasks
            sim.write(open("run%05d.config"%k, "wb"))

        #--- run script
        ips_bin = "runid=$(printf \"%05d\" ${SLURM_NODEID})\n"
        ips_bin += "cd run${runid}\n"
        ips_bin += os.path.join(self.BIN_PATH, self.BIN)
        ips_bin += " --config=../run${runid}.config"
        ips_bin += " --log=../ips_${runid}.log"
        ips_bin += " --platform=../%s"%f_machine_config
        ips_bin += "\ncd .."
        if self.clean_after:
            ips_bin += "\nrm -rf "+rundir
            #ips_bin += "\nrm -f "+rundir+f_machine_config

        print(ips_bin)
        with open("ips_bin.sh", "w") as f: f.write(ips_bin)

        #--- run
        launch_task_option = getattr(self, "LAUNCH_TASK", "")
        if launch_task_option == "srun":
        # This is a temperary implementation to make sure one task per one node on CORI.
            logfile = open("tasks_srun.log", "w")
            errfile = open("tasks_srun.err", "w")
            retcode = subprocess.call(["srun -N %d -n %d -c 64 sh ips_bin.sh"%(ntasks, ntasks)], stdout=logfile, stderr=errfile, shell=True)
        else:
            task_id = services.launch_task(ntasks, cwd, "sh", "ips_bin.sh", logfile="tasks.log", errfile='tasks.err', task_ppn=1)
            retcode = services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing efit')

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        pass
