#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran transport solver component
 -----------------------------------------------------------------------
"""

import os
import shutil
from  component import Component
from Namelist import Namelist
import fastran_io_ps
import fastran_io_instate
import instate_io

class fastran(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timeStamp=0):

        self.icalled = 0
        if not hasattr(self,'USE_INSTATE'): self.USE_INSTATE = 'NO'
        if not hasattr(self,'USE_FASTRAN_NC'): self.USE_FASTRAN_NC = 'NO'
        if not hasattr(self,"PS_BACKEND"): self.PS_BACKEND = "swim"

    def step(self, timeStamp=0):

        #--- entry

        services = self.services

        #--- stage plasma state files

        services.stage_plasma_state()

        #--- get plasma state file name

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        if self.USE_INSTATE == 'YES' or self.PS_BACKEND == 'instate' :
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        if self.USE_FASTRAN_NC == 'YES':
            cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')

        if self.PS_BACKEND == 'swim' :
            cur_bc_file = services.get_config_param('CURRENT_BC')

        #--- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #--- infastran control

        try:
            shutil.copyfile(self.INFASTRAN,"infastran")
        except:
            pass

        infastran = Namelist("infastran")

        for key in ["SOLVE_NE","SOLVE_TE","SOLVE_TI","SOLVE_V","RELAX_J"]:

            iupdate = int(getattr(self,'UPDATE_%s'%key,10000))
            if self.icalled >= iupdate:
                if infastran["infastran"][key][0] == 1:
                   infastran["infastran"][key][0] = 0
                else:
                   infastran["infastran"][key][0] = 1
                print "UPDATE_%s at %d : %d"%(key,self.icalled,infastran["infastran"][key][0])

        # yscale = float(getattr(self,"YSCALE","1.0"))
        # infastran["infastran"]["yscale"] = [yscale]

        infastran.write("infastran")

        #--- generate fastran input

        if self.PS_BACKEND=="instate":
           fastran_io_instate.write_input(cur_instate_file)
        else:
           fastran_io_ps.write_input(cur_state_file, cur_eqdsk_file)

        #--- run fastran

        fastran_bin = os.path.join(self.BIN_PATH, self.BIN)
        print fastran_bin

        ncpu =  int(self.NPROC)
        nky  =  int(self.NPROC_KY)
        n1d  =  ncpu/nky

        print "ncpu = ",ncpu
        print "n1d  = ",n1d
        print "nky  = ",nky

        cwd = services.get_working_dir()
        task_id = services.launch_task(ncpu, cwd, fastran_bin, "%d"%n1d, "%d"%nky, logfile = 'xfastran.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
           print 'Error executing ', 'fastran'
           raise

        #--- update local plasma state

        relax = float(getattr(self,'RELAX',0.5))
        relax_J = float(getattr(self,'RELAX_J',1.0))

        if self.icalled >= int(getattr(self,'ADJUST_IP',10000)):
           adjust_ip = 1
        else:
           adjust_ip = 0

        if self.PS_BACKEND=="instate":
            fastran_io_instate.update_state(cur_instate_file,f_fastran='fastran.nc',relax=relax)
        else:
            fastran_io_ps.update_state(
                f_state=cur_state_file, f_eqdsk=cur_eqdsk_file, f_bc = cur_bc_file,
                f_fastran='fastran.nc', time = timeStamp, relax=relax, relax_J = relax_J, adjust_ip = adjust_ip)

        #--- update plasma state files

        if self.USE_FASTRAN_NC == 'YES':
            shutil.copyfile('fastran.nc',cur_fastran_file)

        if self.USE_INSTATE == 'YES':
            instate_io.ps_to_instate(cur_state_file,cur_eqdsk_file,cur_bc_file,cur_instate_file)

        services.update_plasma_state()

        #--- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        self.icalled = self.icalled + 1
        return

    def finalize(self, timeStamp=0):

        return


def adjust_ip(f_state,f_bc,f_fastran):

    fastran = netCDF4.Dataset(f_fastran,'r',format='NETCDF4')

    ip = fastran.variables["ip"][-1]
    ibs = fastran.variables["ibs"][-1]
    inb = fastran.variables["inb"][-1]
    irf = fastran.variables["irf"][-1]
    fni = (ibs+inb+irf)/ip

    inbc = Namelist(f_bc)
    inbc["inbc"]["ip"][0] = ip/fni
    inbc.write(f_bc)

    print '******* IP ADJUST'
    print ip, ip/fni
