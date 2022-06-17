"""
 -----------------------------------------------------------------------
 massive parallel run
 -----------------------------------------------------------------------
"""
import os
import subprocess
from threading import Thread, Event
from ipsframework.configobj import ConfigObj

from ipsframework import Component
from distributed.diagnostics.plugin import WorkerPlugin


def file_daemon(worker_id, evt, source_dir, target_dir):

    cmd = f"""cd {source_dir}
tar -caf {target_dir}/{worker_id}_archive_$(date +"%Y-%m-%dT%H.%M.%S").tar.gz SUMMARY run?????/*.log run?????/*.out
# Remove all *.tar.gz except most recent 2
ls -t {target_dir}/{worker_id}_archive_*.tar.gz | tail +3 | xargs rm -f
"""

    while not evt.wait(60):  # interval which to archive data
        os.system(cmd)

    os.system(cmd)


class ArchievingPlugin(WorkerPlugin):
    def __init__(self, tmp_dir, target_dir):
        self.tmp_dir = tmp_dir
        self.target_dir = target_dir

    def setup(self, worker):
        self.evt = Event()
        self.thread = Thread(target=file_daemon, args=(worker.id, self.evt, self.tmp_dir, self.target_dir))
        self.thread.start()

    def teardown(self, worker):
        self.evt.set()  # tells the thread to exit
        self.thread.join()


class ips_massive_parallel(Component):
    def step(self, timeid=0):
        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        dir_summary = getattr(self, "SUMMARY", "SUMMARY")
        dir_summary = os.path.realpath(dir_summary)
        if not os.path.exists(dir_summary):
            os.makedirs(dir_summary)

        tmp_xfs = getattr(self, "TMPXFS", "")

        if tmp_xfs:
            tmp_xfs_dir_summary = os.path.realpath(os.path.join(tmp_xfs, "SUMMARY"))

        f_sim_config = getattr(self, "SIMULATION", "")
        sim = ConfigObj(f_sim_config, interpolation='template', file_error=True)

        f_inscan = getattr(self, "INSCAN", "inscan")
        with open(f_inscan,"r") as f:
            inscan = f.readlines()
        nsim = len(inscan)-1
        header = inscan[0]

        dask_nodes = int(getattr(self, "DASK_NODES", "1"))

        try:
            pwd = services.get_config_param("PWD")
        except:
            pwd = os.environ["PWD"]

        #--- add to pool
        pool = services.create_task_pool('pool')
        cwd = services.get_working_dir()

        for k in range(nsim):
            rundir = os.path.realpath(os.path.join(tmp_xfs, "run%05d"%k)) if tmp_xfs else os.path.realpath("run%05d"%k)
            logfile = "ipslog.%05d"%k

            data = inscan[k+1]
            for i, key in enumerate(header.split()):
                print (i, key)
                comp, vname, vtype = key.split(':')
                d = data.split()
                if vtype == 'str':
                    val = str(d[i])
                else:
                    val = eval(vtype+"("+d[i]+")")
                if comp == '':
                    sim[vname] = val
                else:
                    sim[comp][vname] = val

            sim["PWD"] = tmp_xfs if tmp_xfs else pwd
            sim["INPUT_DIR_SIM"] = pwd + "/input"
            sim["RUN_ID"] = "run%05d"%k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(rundir, "run%05d.out"%k)
            sim["USE_PORTAL"] = "False"
            driver = sim['PORTS']['DRIVER']['IMPLEMENTATION']
            sim[driver]["SUMMARY"] = tmp_xfs_dir_summary if tmp_xfs else dir_summary

            services.add_task(
                'pool',
                'task_'+str(k),
                1,
                cwd,
                taskRunner,
                f"run{k:05}",
                str(sim),
                self.TMPXFS,
                logfile=logfile)

        worker_plugin = ArchievingPlugin(self.TMPXFS, dir_summary)

        #--- run
        ret_val = services.submit_tasks('pool', use_dask=True, dask_nodes=dask_nodes,
                                        use_shifter=True,
                                        dask_worker_plugin=worker_plugin)

        print('ret_val = ', ret_val)
        exit_status = services.get_finished_tasks('pool')
        print(exit_status)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)


def taskRunner(runname, sim, tmp_xfs, timeout=1e9):
    if os.system(f'findmnt -nt xfs -T {tmp_xfs}') != 0:
        raise Exception(f"TMPXFS is set but this is either not running in a shifter container "
                        f"or {tmp_xfs} is not mounted as a temporary xfs file")

    rundir = os.path.join(tmp_xfs, runname)

    os.makedirs(rundir, exist_ok=True)
    wrt_localconf(os.path.join(rundir, "local.conf"))

    sim = ConfigObj(eval(sim))
    sim.write(open(os.path.join(rundir, f"{runname}.config"), "wb"))
    
    ips_bin = "which ips.py\n"
    ips_bin += "ips.py"
    ips_bin += f" --config={runname}.config"
    ips_bin += f" --log=ips_{runname}.log"
    ips_bin += " --platform=local.conf"

    ips_bin_path = os.path.join(rundir, "ips_bin.sh")
    with open(ips_bin_path, "w") as f:
        f.write(ips_bin)

    process = subprocess.Popen(["bash", "ips_bin.sh"],
                               cwd=rundir)

    try:
        return process.wait(float(timeout))
    except subprocess.TimeoutExpired:
        process.kill()
        return -1

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