"""
 -----------------------------------------------------------------------
 fastran transport solver component
 -----------------------------------------------------------------------
"""

import os
import shutil
from component import Component
from Namelist import Namelist
import fastran_io_ps
import fastran_io_instate
import instate_io
import zdata
import numpy as np

class fastran(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print ('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('fastran.init() entered')

        self.icalled = 0
        if not hasattr(self, 'USE_INSTATE'): self.USE_INSTATE = 'NO'
        if not hasattr(self, 'USE_FASTRAN_NC'): self.USE_FASTRAN_NC = 'NO'
        if not hasattr(self, 'PS_BACKEND'): self.PS_BACKEND = 'pyps'

    def step(self, timeid=0):
        print('fastran.step() entered')

        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if self.USE_INSTATE == 'YES' or self.PS_BACKEND == 'instate' :
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        if self.USE_FASTRAN_NC == 'YES':
            cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')
        if self.PS_BACKEND == 'pyps' :
            cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- infastran control
        if self.INFASTRAN != "infastran":
            shutil.copyfile(self.INFASTRAN, "infastran")

        infastran = Namelist("infastran")
        for key in ["SOLVE_NE", "SOLVE_TE", "SOLVE_TI", "SOLVE_V", "RELAX_J"]:
            iupdate = int(getattr(self, 'UPDATE_%s'%key, 10000))
            if self.icalled >= iupdate:
                if infastran["infastran"][key][0] == 1:
                   infastran["infastran"][key][0] = 0
                else:
                   infastran["infastran"][key][0] = 1
                print("UPDATE_%s at %d : %d"%(key, self.icalled, infastran["infastran"][key][0]))
        infastran.write("infastran")

        #--- freeze
        ifreeze = int(getattr(self, "FREEZE", 10000))
        irefreeze = int(getattr(self, "REFREEZE", -10000))
        if timeid > ifreeze and timeid < irefreeze:
            print("FASTRAN FREEZE: timeid = %d, ifreeze = %d"%(timeid, ifreeze))
            inmetric_0 = zdata.zdata()
            inmetric_0.read("inmetric")

        #--- generate fastran input
        if self.PS_BACKEND=="instate":
           fastran_io_instate.write_input(cur_instate_file)
        else:
           fastran_io_ps.write_input(cur_state_file, cur_eqdsk_file)

        if timeid > ifreeze and timeid < irefreeze:
            inmetric = zdata.zdata()
            inmetric.read("inmetric")
            inmetric["SHIFT"] = 0.5*np.array(inmetric_0["SHIFT"]) + 0.5*np.array(inmetric["SHIFT"])
            inmetric["PMHD"] = 0.5*np.array(inmetric_0["PMHD"]) + 0.5*np.array(inmetric["PMHD"])
            inmetric.write("inmetric")
        shutil.copyfile('inmetric', 'inmetric_%d'%timeid)

        #--- run fastran
        fastran_bin = os.path.join(self.BIN_PATH, self.BIN)
        print(fastran_bin)

        ncpu = int(self.NPROC)
        nky  = int(self.NPROC_KY)
        n1d  = ncpu/nky

        print("ncpu = ", ncpu)
        print("n1d  = ", n1d)
        print("nky  = ", nky)

        cwd = services.get_working_dir()
        task_id = services.launch_task(ncpu, cwd, fastran_bin, "%d"%n1d, "%d"%nky, logfile = 'xfastran.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           raise Exception('Error executing: fastran')

        #--- update local plasma state
        relax = float(getattr(self, 'RELAX', 0.5))
        relax_J = float(getattr(self, 'RELAX_J', 1.0))
        fni = float(getattr(self, 'FNI', 1.0))
        relax_ip = float(getattr(self, 'RELAX_IP', 0.3))

        if self.icalled >= int(getattr(self, 'ADJUST_IP', 10000)):
           adjust_ip = 1
        else:
           adjust_ip = 0

        if self.PS_BACKEND=="instate":
            fastran_io_instate.update_state(cur_instate_file, f_fastran='fastran.nc', relax=relax)
        else:
            fastran_io_ps.update_state(
                f_state=cur_state_file, f_eqdsk=cur_eqdsk_file, f_bc=cur_bc_file,
                f_fastran='fastran.nc', time = timeid, relax=relax, relax_J=relax_J, adjust_ip=adjust_ip, fni_target=fni, relax_ip=relax_ip)

        #--- update plasma state files
        if self.USE_FASTRAN_NC == 'YES':
            shutil.copyfile('fastran.nc', cur_fastran_file)

        if self.USE_INSTATE == 'YES':
            instate_io.ps_to_instate(cur_state_file, cur_eqdsk_file, cur_bc_file, cur_instate_file)

        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

        self.icalled = self.icalled + 1

    def finalize(self, timeid=0):
        print('fastran.finalize() entered')

        services = self.services

        cwd = services.get_working_dir()
        ishot = int(services.get_config_param('SHOT_NUMBER'))
        itime = int(services.get_config_param('TIME_ID'))

        dir_state = services.get_config_param('PLASMA_STATE_WORK_DIR')
        f = "f%06d.%05d"%(ishot, itime)
        shutil.copyfile("fastran.nc", os.path.join(dir_state, f))

def adjust_ip(f_state, f_bc, f_fastran):

    fastran = netCDF4.Dataset(f_fastran,'r',format='NETCDF4')

    ip = fastran.variables["ip"][-1]
    ibs = fastran.variables["ibs"][-1]
    inb = fastran.variables["inb"][-1]
    irf = fastran.variables["irf"][-1]
    fni = (ibs+inb+irf)/ip

    inbc = Namelist(f_bc)
    inbc["inbc"]["ip"][0] = ip/fni
    inbc.write(f_bc)

    print ('******* IP ADJUST')
    print (ip, ip/fni)
