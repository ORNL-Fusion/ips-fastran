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

from ipsframework import Component
from ipsframework import ipsutil


class ips_scan(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print("Created %s" % (self.__class__))

    def init(self, timeid=0):
        self.clean_after = int(getattr(self, "CLEAN_AFTER", "0"))
        self.time_out = int(getattr(self, "TIME_OUT", "3600000"))

    def step(self, timeid=0):
        # --- stage plasma state files
        self.services.stage_state()

        # --- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # --- setup
        dir_summary = getattr(self, "SUMMARY", "")
        if dir_summary:
            dir_summary = os.path.realpath(dir_summary)
            if not os.path.exists(dir_summary):
                os.makedirs(dir_summary)

        f_sim_config = getattr(self, "SIMULATION", "")
        sim = ConfigObj(f_sim_config, interpolation="template", file_error=True)

        f_inscan = getattr(self, "INSCAN", "inscan")
        with open(f_inscan, "r") as f:
            inscan = f.readlines()
        nsim = len(inscan) - 1
        header = inscan[0]

        use_dask = bool(int(getattr(self, "USE_DASK", "0")))
        dask_nodes = int(getattr(self, "DASK_NODES", "1"))
        task_ppn = int(getattr(self, "TASK_PPN", ""))

        rank = int(getattr(self, "RANK", "0"))
        size = int(getattr(self, "SIZE", "1"))
        print(f"NSIM = {nsim} at RANK = {rank}")

        try:
            pwd = self.services.get_config_param("PWD")
        except BaseException:
            pwd = os.environ["PWD"]

        scan_dir = getattr(self, "SCAN_DIR", "")
        if scan_dir:
            if not os.path.exists(scan_dir):
                os.makedirs(scan_dir)
        base_dir = scan_dir if scan_dir else cwd

        cwd = self.services.get_working_dir()
        wrt_localconf(os.path.join(base_dir, "local.conf"))

        # --- add to pool
        pool = self.services.create_task_pool("pool")
        tasks = {}
        for k in range(rank, nsim, size):
            if scan_dir:
                rundir = os.path.realpath(os.path.join(scan_dir, "run%05d" % k))
                logfile = os.path.realpath(os.path.join(scan_dir, "log%05d" % k))
            else:
                rundir = os.path.realpath("run%05d" % k)
                logfile = "log%05d" % k
            if not os.path.exists(rundir):
                os.makedirs(rundir)

            data = inscan[k + 1]
            for i, key in enumerate(header.split()):
                print(i, key)
                comp, vname, vtype = key.split(":")
                d = data.split()
                if vtype == "str":
                    val = str(d[i])
                else:
                    val = eval(vtype + "(" + d[i] + ")")
                if comp == "":
                    sim[vname] = val
                else:
                    sim[comp][vname] = val

            sim["PWD"] = pwd
            sim["RUN_ID"] = "run%05d" % k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(base_dir, "run%05d.out" % k)
            sim["LOG_FILE"] = os.path.join(base_dir, "run%05d.log" % k)
            try:
                sim["PARENT_PORTAL_RUNID"] = self.services.get_config_param("PORTAL_RUNID")
                sim["USE_PORTAL"] = "True"
            except Exception:
                sim["USE_PORTAL"] = "False"
            driver = sim["PORTS"]["DRIVER"]["IMPLEMENTATION"]
            # sim[driver]["SUMMARY"] = dir_summary
            sim.write(open(os.path.join(base_dir, "run%05d.config" % k), "wb"))

            dir_state = os.path.join(rundir, "work/plasma_state")  # not used
            dir_summary_local = os.path.join(base_dir, f"summary{k:05}")
            f_tarfile = os.path.join(dir_summary, f"summary{k:05}.tar")

            ips_bin = os.path.join(self.BIN_PATH, self.BIN)
            ips_bin += " --config=run%05d.config" % k
            ips_bin += " --log=ips%05d.log" % k
            ips_bin += " --platform=local.conf"
            ips_bin += f"\ncollect.py --rdir0={rundir} --rdir={rundir} --sdir={dir_summary_local} --time={k} --single --input={os.path.join(cwd,'collect.json')}"
            # ips_bin += f"\ncp {f_eventlog} {f_eventlog_target}"
            ips_bin += f"\ncp {os.path.join(rundir, 'simulation_log/*.json')} {os.path.join(dir_summary_local, f'eventlog_{k:05}')}"
            ips_bin += f"\ncp run{k:05}.log {dir_summary_local}/"
            ips_bin += f"\ncp run{k:05}.out {dir_summary_local}/"
            ips_bin += f"\ntar -C {dir_summary_local} -cvf {f_tarfile} ."
            if self.clean_after:
                ips_bin += "\nrm -rf " + rundir

            with open(os.path.join(rundir, "ips_bin.sh"), "w") as f:
                f.write(ips_bin)

            self.services.add_task(
                "pool",
                "task_" + str(k),
                1,
                base_dir,
                "sh",
                os.path.join(rundir, "ips_bin.sh"),
                logfile=logfile,
                timeout=self.time_out,
                task_ppn=task_ppn)

        # --- run
        ret_val = self.services.submit_tasks("pool", use_dask=use_dask, dask_nodes=dask_nodes)
        print("ret_val = ", ret_val)
        exit_status = self.services.get_finished_tasks("pool")
        print(exit_status)

        # --- update plasma state files
        self.services.update_state()

        # --- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
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
    f = open(fname, "w")
    f.write(s)
    f.close()
