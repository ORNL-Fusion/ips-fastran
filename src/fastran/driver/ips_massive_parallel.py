"""
 -----------------------------------------------------------------------
 massive parallel run
 -----------------------------------------------------------------------
"""
import os
import shutil
import subprocess
from threading import Thread, Event
from configobj import ConfigObj

from ipsframework import Component
from distributed.diagnostics.plugin import WorkerPlugin


def file_daemon(worker_id, evt, source_dir, target_dir, archive_files):
    files = archive_files if archive_files else "run?????/run?????.???  run?????/simulation_log/*.eventlog"
    files = "SUMMARY " + files
    print(">>> archive_files:", files)

    cmd = f"""cd {source_dir}
tar -caf {target_dir}/{worker_id}_archive_$(date +"%Y-%m-%dT%H.%M.%S").tar.gz {files}
# Remove all *.tar.gz except most recent 2
ls -t {target_dir}/{worker_id}_archive_*.tar.gz | tail +3 | xargs rm -f
"""
    print(">>> cmd:", cmd)

#        cmd = f"""cd {source_dir}
#tar -caf {target_dir}/{worker_id}_archive_$(date +"%Y-%m-%dT%H.%M.%S").tar.gz SUMMARY run?????/run?????.???  run?????/simulation_log/*.eventlog
## Remove all *.tar.gz except most recent 2
#ls -t {target_dir}/{worker_id}_archive_*.tar.gz | tail +3 | xargs rm -f
#"""

    while not evt.wait(600):  # interval which to archive data
        os.system(cmd)

    os.system(cmd)


class ArchievingPlugin(WorkerPlugin):
    def __init__(self, tmp_dir, target_dir, archive_files):
        self.tmp_dir = tmp_dir
        self.target_dir = target_dir
        self.archive_files = archive_files

    def setup(self, worker):
        self.evt = Event()
        self.thread = Thread(target=file_daemon, args=(worker.id, self.evt, self.tmp_dir, self.target_dir, self.archive_files))
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

        task_ppn = int(getattr(self, "TASK_PPN", "-1"))
        task_ppn = task_ppn if task_ppn > 0 else services.get_config_param("PROCS_PER_NODE")
        print('task_ppn =', task_ppn)

        task_nproc = int(getattr(self, "TASK_NPROC", "1"))

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
           #sim["INPUT_DIR_SIM"] = pwd + "/input"
            sim["INPUT_DIR_SIM_TEMP"] = os.path.realpath(os.path.join(pwd , "input"))
            sim["INPUT_DIR_SIM"] = os.path.join(rundir , "input")
            sim["RUN_ID"] = "run%05d"%k
            sim["SIM_ROOT"] = rundir
            sim["OUT_REDIRECT"] = "True"
            sim["OUT_REDIRECT_FNAME"] = os.path.join(rundir, "run%05d.out"%k)
            sim["LOG_FILE"] = os.path.join(rundir, "run%05d.log"%k)
            sim["USE_PORTAL"] = "True"
            try:
                sim['PARENT_PORTAL_RUNID'] = self.services.get_config_param("PORTAL_RUNID")
                sim["USE_PORTAL"] = "True"
            except Exception:
                sim["USE_PORTAL"] = "False"
            sim["TASK_NPROC"] = task_nproc
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

        archive_files = getattr(self, "ARCHIVE_FILES", "")
        worker_plugin = ArchievingPlugin(self.TMPXFS, dir_summary, archive_files)

        #--- run
        ret_val = services.submit_tasks('pool', use_dask=True, dask_nodes=dask_nodes,
                                        use_shifter=True,
                                        dask_ppw=task_ppn,
                                        dask_worker_plugin=worker_plugin)

        print('ret_val = ', ret_val)
        exit_status = services.get_finished_tasks('pool')
        services.remove_task_pool('pool')
        print(exit_status)

        clean_after = int(getattr(self, 'CLEAN_AFTER', '0'))
        if clean_after:
            with open('cmd.sh', 'w') as f:
                 f.write(f'rm -rf {tmp_xfs}/run?????\n')

            cwd = self.services.get_working_dir()
            cmd =  f'shifter sh cmd.sh'
            task_id = self.services.launch_task(self.DASK_NODES, cwd, cmd, task_ppn=1, logfile='clean.log')
            retcode = self.services.wait_task(task_id)
            if (retcode != 0):
                e = 'Error executing command:  clean '
                raise Exception(e)

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

    sim = ConfigObj(eval(sim))

    task_nproc = int(sim["TASK_NPROC"])

    input_dir = sim["INPUT_DIR_SIM_TEMP"]
    print(input_dir )
    shutil.copytree(input_dir, os.path.join(rundir, "input"))
    #sim["INPUT_DIR_SIM"] =  os.path.join(rundir, "input")

    sim.write(open(os.path.join(rundir, f"{runname}.config"), "wb"))

#   ips_bin = "which ips.py\n"
    ips_bin = "echo $PYTHONPATH\n"
    ips_bin += "which ips.py\n"
    ips_bin += "ips.py"
    ips_bin += f" --config={runname}.config"
    ips_bin += f" --log=ips_{runname}.log"
    ips_bin += " --platform=local.conf"

    ips_bin_path = os.path.join(rundir, "ips_bin.sh")
    with open(ips_bin_path, "w") as f:
        f.write(ips_bin)

    wrt_localconf(os.path.join(rundir, "local.conf"), task_nproc)

    process = subprocess.Popen(["bash", "ips_bin.sh"],
                               cwd=rundir)

    try:
        return process.wait(float(timeout))
    except subprocess.TimeoutExpired:
        process.kill()
        return -1

def wrt_localconf(fname="local.conf", task_nproc=1):
    if task_nproc == 1:
       cmd = "eval"
    else:
       cmd = "mpiexec"
    s = \
"""HOST = local
MPIRUN = %s
NODE_DETECTION = manual
PROCS_PER_NODE = %d
CORES_PER_NODE = %d
SOCKETS_PER_NODE = 1
NODE_ALLOCATION_MODE = SHARED
USE_ACCURATE_NODES = ON
"""%(cmd, task_nproc, task_nproc)

    f=open(fname, "w")
    f.write(s)
    f.close()
