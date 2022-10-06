"""
 -----------------------------------------------------------------------
 fastran transport solver component
 -----------------------------------------------------------------------
"""

import os
import shutil
import numpy as np
from ipsframework import Component
from Namelist import Namelist
from fastran.solver import fastran_io_ps
from fastran.solver import fastran_io_instate
from fastran.instate import instate_io
from fastran.solver import zdata


class fastran(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('fastran.init() entered')
        self.icalled = 0

    def step(self, timeid=0):
        print('fastran.step() entered')
        ifreeze = int(getattr(self, "FREEZE", -1))
        iresume = int(getattr(self, "RESUME", -1))
        if ifreeze >= 0 and timeid >= ifreeze:
            if iresume < 0 or timeid < iresume:
                print("fastran skipped, FREEZE = %d, RESUME = %d, TIMEID = %d" % (ifreeze, iresume, timeid))
                return None

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file name
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        ps_backend = getattr(self, 'PS_BACKEND', 'pyps').lower()
        update_instate = getattr(self, 'UPDATE_INSTATE', 'enabled').lower()

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- infastran control
        if self.INFASTRAN != "infastran":
            shutil.copyfile(self.INFASTRAN, "infastran")

        infastran = Namelist("infastran")
        for key in ["SOLVE_NE", "SOLVE_TE", "SOLVE_TI", "SOLVE_V", "RELAX_J"]:
            iupdate = int(getattr(self, 'UPDATE_%s' % key, -1))
            if self.icalled >= iupdate and iupdate >= 0:
                if infastran["infastran"][key][0] == 1:
                    infastran["infastran"][key][0] = 0
                else:
                    infastran["infastran"][key][0] = 1
                print("UPDATE_%s at %d : %d" % (key, self.icalled, infastran["infastran"][key][0]))
        infastran.write("infastran")

        # -- generate fastran input
        if ps_backend == "instate":
            fastran_io_instate.write_input(cur_instate_file)
        else:
            recycle = float(getattr(self, "RECYCLE", "0"))
            print("recycle =", recycle)
            fastran_io_ps.write_input(cur_state_file, cur_eqdsk_file, recycle=recycle)

        # -- run fastran
        fastran_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(fastran_bin)

        ncpu = int(self.NPROC)
        nky = int(self.NPROC_KY)
        n1d = ncpu//nky

        print("ncpu = ", ncpu)
        print("n1d  = ", n1d)
        print("nky  = ", nky)

        cwd = self.services.get_working_dir()
        task_id = self.services.launch_task(ncpu, cwd, fastran_bin, "%d" % n1d, "%d" % nky, logfile='xfastran.log')
        retcode = self.services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing: fastran')

        # -- update local plasma state
        relax = float(getattr(self, 'RELAX', 0.5))
        relax_J = float(getattr(self, 'RELAX_J', 1.0))
        fni = float(getattr(self, 'FNI', 1.0))
        relax_ip = float(getattr(self, 'RELAX_IP', 0.3))

        if self.icalled >= int(getattr(self, 'ADJUST_IP', 10000)):
            adjust_ip = 1
        else:
            adjust_ip = 0

        if ps_backend == "instate":
            fastran_io_instate.update_state(
                f_instate=cur_instate_file, 
                f_fastran='fastran.nc', 
                relax=relax)
        else:
            fastran_io_ps.update_state(
                f_state=cur_state_file, 
                f_eqdsk=cur_eqdsk_file, 
                f_instate=cur_instate_file,
                f_fastran='fastran.nc', 
                time=timeid, 
                relax=relax, 
                relax_J=relax_J, 
                adjust_ip=adjust_ip, 
                fni_target=fni, 
                relax_ip=relax_ip)

        if update_instate == 'enabled' and ps_backend == 'pyps':
            print("updating instate")
            instate_io.ps_to_instate(cur_state_file, cur_eqdsk_file, cur_instate_file)

        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

        self.icalled = self.icalled + 1

    def finalize(self, timeid=0):
        print('fastran.finalize() entered')

        cwd = self.services.get_working_dir()
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(self.services.get_config_param('TIME_ID'))

        dir_state = self.services.get_config_param('STATE_WORK_DIR')
        f = "f%06d.%05d" % (ishot, itime)
        shutil.copyfile("fastran.nc", os.path.join(dir_state, f))


def adjust_ip(f_state, f_bc, f_fastran):

    fastran = netCDF4.Dataset(f_fastran, 'r', format='NETCDF4')

    ip = fastran.variables["ip"][-1]
    ibs = fastran.variables["ibs"][-1]
    inb = fastran.variables["inb"][-1]
    irf = fastran.variables["irf"][-1]
    fni = (ibs+inb+irf)/ip

    inbc = Namelist(f_bc)
    inbc["inbc"]["ip"][0] = ip/fni
    inbc.write(f_bc)

    print('******* IP ADJUST')
    print(ip, ip/fni)
