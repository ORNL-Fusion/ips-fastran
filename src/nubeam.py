#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 nubeam component for steady-state solution
 JM
 -----------------------------------------------------------------------
"""

import sys,os,os.path,shutil,pickle,glob
import subprocess
from numpy import *
import netCDF4

from  component import Component

# import zcode libraries
import Namelist
from zplasmastate import plasma_state_file
import znubeam
from plasmastate import *
import traceback

class nubeam(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def set_env_vars(self):

        services = self.services
        workdir = services.get_working_dir()

        try:
            os.environ['ADASDIR'] = self.ADAS
        except Exception:
            services.exeception('no ADAS parameter')
            raise

        try:
            preact =  self.PREACT
            print "preact = ", preact
        except Exception:
            services.exeception('no PREACT parameter')
            raise

        try:
            copy_preact = self.COPY_PREACT
        except:
            copy_preact = 0

        if copy_preact:
            try:
                print 'Copying PREACT dir ...\n'
                print 'Source : ' + preact + '\n'
                print 'Destination : ' + os.path.join(workdir, "PREACT") + '\n'
                shutil.copytree(preact, os.path.join(workdir, "PREACT"))
            except:
                print 'PREACT directory not copied or already there.\n'
            os.environ['PREACTDIR'] = 'PREACT'
        else:
            print 'Setting PREACTDIR env var to ' + self.PREACT + '\n'    
            os.environ['PREACTDIR'] = self.PREACT

        print 'PREACTDIR : ' + os.environ['PREACTDIR'] + '\n'    

        return

    def init(self, timeStamp=0.0):

        #--- entry

        services = self.services

        services.stage_plasma_state()

        pstool_path = self.PSTOOL_PATH

        pstool_bin = pstool_path

        print 'pstool_bin path:'
        print pstool_bin

        print 'pstool nubeam init'

        logfile=open("pstool_init.log","w")
        retcode = subprocess.call([pstool_bin, "init", "nubeam"],
                      stdout=logfile,stderr=logfile)
        logfile.close()
        if (retcode != 0):
            logMsg = 'Error executing pstool_bin'
            services.exception(logMsg)
            raise Exception(logMgs, pstool_bin)

        cur_state_file = services.get_config_param('CURRENT_STATE')
        shutil.copyfile("ips-state-nubeam.nc",cur_state_file)

        try:
            services.update_plasma_state()
        except Exception:
            services.exception('Error in call to update_plasma_state()')
            raise

        #--- get plasma state file name

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)
        
        #--- set environment vars

        self.set_env_vars()

        #--- generate nubeam input
 
        nubeam_files = Namelist.Namelist()
        nubeam_files["NUBEAM_FILES"]["INPUT_PLASMA_STATE"] = \
            [cur_state_file]
        #nubeam_files["NUBEAM_FILES"]["OUTPUT_PLASMA_STATE"] = \
        #    [cur_state_file]
        nubeam_files["NUBEAM_FILES"]["PLASMA_STATE_UPDATE"] = \
            ["state_changes.cdf"]
        nubeam_files["NUBEAM_FILES"]["INIT_NAMELIST"] = \
            ["nubeam_init_input.dat"]
        nubeam_files.write("nubeam_init_files.dat")

        nubeam_files = Namelist.Namelist()
        nubeam_files["NUBEAM_FILES"]["INPUT_PLASMA_STATE"] = \
            [cur_state_file]
        #nubeam_files["NUBEAM_FILES"]["OUTPUT_PLASMA_STATE"] = \
        #    [cur_state_file]
        nubeam_files["NUBEAM_FILES"]["PLASMA_STATE_UPDATE"] = \
            ["state_changes.cdf"]
        nubeam_files["NUBEAM_FILES"]["STEP_NAMELIST"] = \
            ["nubeam_step_input.dat"]
        nubeam_files.write("nubeam_step_files.dat")

        innubeam = Namelist.Namelist("innubeam")

        nubeam_init_input = Namelist.Namelist()
        nubeam_init_input["NBI_INIT"] = innubeam["NBI_INIT"]
        nubeam_init_input.write("nubeam_init_input.dat")

        nubeam_step_input = Namelist.Namelist()
        nubeam_step_input["NBI_UPDATE"] = innubeam["NBI_UPDATE"]
        nubeam_step_input.write("nubeam_step_input.dat")

        #--- stage plasma state files

        services.stage_plasma_state()

        try:
            shutil.copyfile(cur_state_file, "input_state.cdf")
        except Exception:
            logMsg = 'Error copying plasma state files over to generic file names'     
            services.exception(logMsg)
            raise 

        #--- setup nubeam_comp_exec run

        os.environ['NUBEAM_ACTION'] = 'init'
        os.environ['FRANTIC_INIT'] = '50'
        try:
            del os.environ['FRANTIC_ACTION']
        except:
            pass

        try:
            nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            nubeam_bin = os.path.join(self.BIN_PATH, 'mpi_nubeam_comp_exec')

        print "NUBEAM Binary: "
        print nubeam_bin

        task_id = services.launch_task(1, workdir, nubeam_bin,
                      logfile = 'log.init_nubeam')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            logMsg = 'Error executing command:  mpi_nubeam_comp_exec: init '
            services.exception(logMsg)
            raise Exception(logMsg)

        try:
            shutil.copyfile('log.init_nubeam', 'log.nubeam')
        except Exception:
            logMsg = 'Error in file cp'
            services.exception(logMsg)
            raise 

        #--- update plasma state

        # services.update_plasma_state() #<== ugly
        services.merge_current_plasma_state("state_changes.cdf",
                     logfile='log.update_state')

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return

    def step(self, timeStamp=0.0):

        if (self.services == None) :
            logMsg = 'Error in nubeam.step() : No services'
            raise Exception(logMsg)

        services = self.services

        #--- get plasma state file name

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage plasma state files

        try:
            services.stage_plasma_state()
            shutil.copyfile (cur_state_file,'initial-'+cur_state_file)
        except Exception:
            logMsg = 'Error in call to stage_plasma_state()'
            services.exception(logMsg)
            raise

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- get work directory

        workdir = services.get_working_dir()

        #--- set nubeam_comp_exec run

        innubeam = Namelist.Namelist("innubeam")

        # Santiy check on plasma state density

        _ps = PlasmaState("ips",1)
        _ps.read(cur_state_file)
        nS = _ps["ns"].shape[0]

        for nn in range(0,nS-1):

            print 'Checking for NaNs before running Nubeam'

            if np.isnan(np.sum(_ps["ns"][nn])):

                print 'ERROR : NaN detected before running Nubeam'
                traceback.print_stack()
                services.error('ERROR : NaN detected before running Nubeam')
                raise

        ncpu =  int(self.NPROC)
        nstep = innubeam["nubeam_run"]["nstep"][0]
        try:
            navg = innubeam["nubeam_run"]["navg"][0]
        except:
            navg = 0
        if nstep-navg < 0: 
           raise Exception("nubeam.py: nstep < navg")

        dt_nubeam = innubeam["nubeam_run"]["dt_nubeam"][0]

        print "ncpu  = ",ncpu
        print "dt    = ",dt_nubeam
        print "nstep = ",nstep
        print "navg  = ",navg

        difb_0  = innubeam["nbi_model"]["difb_0"  ][0]
        difb_a  = innubeam["nbi_model"]["difb_a"  ][0]
        difb_in = innubeam["nbi_model"]["difb_in" ][0]
        difb_out= innubeam["nbi_model"]["difb_out"][0]

        ps = plasma_state_file(cur_state_file)
        rho_anom = ps["rho_anom"][:]
        ps["difb_nbi"][:] = difb_a \
            +(difb_0-difb_a)*(1.0-rho_anom**difb_in)**difb_out
        ps.close()

        services.update_plasma_state()

        try:
            nubeam_bin = os.path.join(self.BIN_PATH, self.BIN)
        except:
            nubeam_bin = os.path.join(self.BIN_PATH, 'mpi_nubeam_comp_exec')

        os.environ['NUBEAM_ACTION'] = 'STEP'
        os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f'%(nstep-navg,dt_nubeam)
        os.environ['STEPFLAG'] = 'TRUE'
        os.environ['NUBEAM_POSTPROC'] = 'summary_test' #'FBM_WRITE'
        try:
            del os.environ['FRANTIC_INIT']
        except:
            pass
        os.environ['FRANTIC_ACTION'] = 'execute' #'none' #'execute'

        #--- run nstep-navg

        print os.environ['NUBEAM_REPEAT_COUNT'] 
        task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, logfile = 'log.nubeam')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
            print nubeam_bin
            logMsg = 'Error executing command:  mpi_nubeam_comp_exec: step '
            raise Exception(logMsg)

        if navg > 0:

            print 'run navg'

             #--- run navg

            os.environ['NUBEAM_REPEAT_COUNT'] = '%dx%f'%(1,dt_nubeam)
            for k in range(navg):
                task_id = services.launch_task(self.NPROC, workdir, nubeam_bin, logfile = 'log.nubeam2')
                retcode = services.wait_task(task_id)
                if (retcode != 0):
                    logMsg = 'Error executing command:  mpi_nubeam_comp_exec: step avg '
                    raise Exception(logMsg)
                shutil.copyfile("state_changes.cdf","state_changes_%d.cdf"%k)

             #--- avgerage

            data = {}
            for k in range(navg):
                filename = "state_changes_%d.cdf"%k
                data[k] = netCDF4.Dataset(filename,'r+',format='NETCDF4')

            vars = []
            for v in data[0].variables.keys():
                if v not in ["ps_partial_update", "version_id"]: vars.append(v)

            for var in vars:
                avg = []
                for k in range(navg):
                    avg.append(data[k].variables[var][:])
                avg = average(array(avg),axis=0)
                data[0].variables[var][:] =  avg
            
            for k in range(navg):
                data[k].close()

            shutil.copyfile ("state_changes_0.cdf","state_changes.cdf")

        #--- update plasma state

        # services.update_plasma_state() 

        services.merge_current_plasma_state("state_changes.cdf", logfile='log.update_state')

        services.stage_plasma_state()
        znubeam.update_ps_profile(cur_state_file,cur_eqdsk_file)
        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        return


    def restart(self, timeStamp=0.0):

        services = self.services
        componentName = self.__class__.__name__ 
        logMsg = 'INFO : restart() called for '+componentName
        services.info(logMsg)
        restart_root = services.get_config_param('RESTART_ROOT')
        services.get_restart_files(restart_root, timeStamp, self.RESTART_FILES)

        #--- set environment vars

        self.set_env_vars()

        return
    

    def checkpoint(self, timestamp=0.0):

        services = self.services

        services.info('checkpoint() called for nubeam')

        services.save_restart_files(timestamp, self.RESTART_FILES)

        return

    def finalize(self, timeStamp=0.0):

        return
