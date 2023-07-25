"""
 -----------------------------------------------------------------------
 chease component
 -----------------------------------------------------------------------
"""

import os
import re
import sys
import glob
import time
import shutil
import subprocess

from time         import time
from numpy        import mean,pi,sqrt,array,linspace,argmin,argmax
from Namelist     import Namelist
from datetime     import datetime
from ipsframework import Component

from scipy.interpolate   import interp1d
from fastran.equilibrium import efit_io
from fastran.equilibrium import efit_eqdsk
from fastran.equilibrium import cheasefiles
from fastran.equilibrium import cheasetools
from fastran.equilibrium import cheasewrapper

from fastran.instate                 import instate_io
from fastran.solver.inmetric_io      import ps_to_inmetric
from fastran.equilibrium.efit_eqdsk  import readg, writeg, scaleg
from fastran.plasmastate.plasmastate import plasmastate

class cheasepy(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('cheasepy.init() started')
        services = self.services

        #--- get shot and time
        self.TOKAMAK = services.get_config_param('TOKAMAK_ID')
        self.SHOT_ID = int(services.get_config_param('SHOT_NUMBER'))
        self.TIME_ID = int(self.TIME_ID if hasattr(self, "TIME_ID") else services.get_config_param('TIME_ID'))

        #--- get plasma state file names
        self.plasma_state_dir = services.get_config_param('STATE_WORK_DIR')

        self.boundfname       = services.get_config_param('CURRENT_BC')
        self.eqdskfname       = services.get_config_param('CURRENT_EQDSK')
        self.statefname       = services.get_config_param('CURRENT_STATE')
        self.instatefname     = services.get_config_param('CURRENT_INSTATE')

        #--- get path to inputs
        self.inputfpath = getattr(self,'INPUT_DIR',  '')

        print('cheasepy.init() done')

    def step(self, timeid=0):
        print('cheasepy.step() started')
        #--- code entry
        services = self.services

        #--- GENE EXCUTABLE
        chease_bin = os.path.join(self.BIN_PATH, self.BIN)
        print("CHEASE Binary: %s" % chease_bin)

        #--- generate inefit
        mode = getattr(self, "PRESSURE", "kinetic")

        setParam = {}
        setParam['mode'] = mode

        #--- get working directory
        self.cwd = services.get_working_dir()

        #--- stage plasma state files
        services.stage_state()

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        ps_backend  =     getattr(self, 'PS_BACKEND',  'PS')
        init_run    = int(getattr(self, 'INIT_RUN',    "0"))
        eq_recycle  = int(getattr(self, 'EQ_RECYCLE',  "0"))

        setParam['init_run'] = init_run

        #--- read input files
        self.inputfname = getattr(self,'INPUT_FILES','inchease').split()
        if   'inchease0' in self.inputfname:
              self.incheasefname = 'inchease0'
        elif 'inchease' in self.inputfname:
            self.incheasefname = 'inchease'
        else:
            print("IOError: inchease file is missing. EXIT!"); sys.exit()

        self.path_to_inchease = os.path.join(self.inputfpath,self.incheasefname)
        if os.path.isfile(self.path_to_inchease):
           inchease = Namelist(self.path_to_inchease, case="lower")
        else:
           inchease = {}

        namelistVals = {}
        if 'namelist' in inchease.keys():
            for ikey in inchease['namelist'].keys():
                if ikey in ['cocos_in','cocos_out']:
                    namelistVals[ikey] = inchease['namelist'][ikey][0]
                else:
                    namelistVals[ikey.upper()] = inchease['namelist'][ikey][0]

        if   'NRBOX'    not in namelistVals: namelistVals['NRBOX']    = 129
        if   'NZBOX'    not in namelistVals: namelistVals['NZBOX']    = 129
        if   'RELAX'    not in namelistVals: namelistVals['RELAX']    = 0.0
        if   'NFUNRHO'  not in namelistVals: namelistVals["NFUNRHO"]  = 1
        if   'NRHOMESH' not in namelistVals: namelistVals["NRHOMESH"] = 1

        importedVals = {}
       #if 'imported'  in inchease.keys():
       #    for ikey in inchease['imported'].keys():
       #        importedVals[ikey] = inchease['imported'][ikey][0]

        srcVals = {}
        if 'sources'  in inchease.keys():
            for ikey in inchease['sources'].keys():
                srcVals[ikey] = inchease['sources'][ikey][0]

        #--- set default for source values
        if 'rhomesh_src'   not in srcVals: srcVals['rhomesh_src']   = 'imported'
        if 'current_src'   not in srcVals: srcVals['current_src']   = 'imported'
        if 'pressure_src'  not in srcVals: srcVals['pressure_src']  = 'imported'
        if 'boundary_src'  not in srcVals: srcVals['boundary_src']  = 'imported'
        if 'eprofiles_src' not in srcVals: srcVals['eprofiles_src'] = 'imported'
        if 'iprofiles_src' not in srcVals: srcVals['iprofiles_src'] = 'imported'
        if 'boundary_type' not in srcVals: srcVals['boundary_type'] = 'asis'
       #if init_run:
       #    if 'boundary_type' not in srcVals: srcVals['boundary_type'] = 'interp'
       #else:
       #    if 'boundary_type' not in srcVals: srcVals['boundary_type'] = 'asis'

        if 'instate' in self.inputfname:
            self.instatefname = 'instate'
            self.path_to_instate = os.path.join(self.inputfpath,self.instatefname)

        #--- create input plasma profiles
        if   ps_backend == "PS":
            if init_run:
             statedata = cheasefiles.get_plasmastate(statefpath=self.statefname,bcfpath=self.boundfname,setParam=setParam)
            else:
             statedata = cheasefiles.get_plasmastate(statefpath=self.statefname,eqfpath=self.eqdskfname,setParam=setParam)
        elif ps_backend == "INSTATE":
             statedata = cheasefiles.get_plasmastate(instatefpath=self.instatefname,bcfpath=self.boundfname,setParam=setParam)

        mu0 = 4.0e-7*pi

        importedVals['ITEXP']     = statedata['ip']
        importedVals['R0EXP']     = statedata['RCTR']
        importedVals['B0EXP']     = statedata['BCTR']
        importedVals['rbound']    = statedata['rbound']
        importedVals['zbound']    = statedata['zbound']

        namelistVals['NIDEAL']    = 6
        namelistVals['R0EXP']     = statedata['RCTR']
        namelistVals['B0EXP']     = statedata['BCTR']
        namelistVals['CURRT']     = statedata['ip']*mu0/statedata['RCTR']/statedata['BCTR'] 
        namelistVals['PSIBNDEXP'] = 0.0

        importedVals['Te']        = statedata['Te']
        importedVals['ne']        = statedata['ne']
        importedVals['Ti']        = statedata['Ti']
        importedVals['ni']        = statedata['ni']
        importedVals['nz']        = statedata['nz']
        importedVals['Zeff']      = statedata['zeff']
        importedVals['rhotor']    = statedata['rho']

        if init_run == 1:
            importedVals['ZMAX']     = (max(statedata['zbound']) + min(statedata['zbound'])) / 2.0
           #importedVals['ZMAX']     = 0.0 # (max(statedata['zbound']) + min(statedata['zbound'])) / 2.0

           #namelistVals['RC']       = 1.0
           #namelistVals['R0']       = 1.0
           #namelistVals['RZ0']      = 0.0
            namelistVals['ZBOXMID']  = 0.0
            namelistVals['ZBOXLEN']  = statedata['ZLEN']
            namelistVals['RBOXLFT']  = statedata['RLFT']
            namelistVals['RBOXLEN']  = statedata['RLEN']

            namelistVals['NCSCAL']   = 1
            namelistVals['PREDGE']   = statedata['pressure'][-1]

            importedVals['Jprl']     = statedata['jpar']
            importedVals['pressure'] = statedata['pressure']

        else:
          # eqdskdata = cheasefiles.read_eqdsk(fpath=os.path.join(self.plasma_state_dir,self.eqdskfname))
          # importedVals['ZMAX']     = eqdskdata['ZMAX'] # 0.0 # (max(statedata['zbound']) + min(statedata['zbound'])) / 2.0
          # namelistVals['ZBOXMID']  = eqdskdata['ZMID']    # 0.0
          # namelistVals['ZBOXLEN']  = eqdskdata['ZLEN']
          # namelistVals['RBOXLFT']  = eqdskdata['RLFT']
          # namelistVals['RBOXLEN']  = eqdskdata['RLEN']

          # namelistVals['NCSCAL']   = 4
          # namelistVals['PREDGE']   = cheasetools.interp(eqdskdata['rhopsi'],eqdskdata['pressure'],statedata['rho'])[-1]

          # importedVals['Jprl']     = cheasetools.interp(eqdskdata['rhopsi'],eqdskdata['jtot'],statedata['rho'])
          # importedVals['pressure'] = cheasetools.interp(eqdskdata['rhopsi'],eqdskdata['pressure'],statedata['rho'])

            importedVals['ZMAX']     = 0.0 # (max(statedata['zbound']) + min(statedata['zbound'])) / 2.0
           #namelistVals['RC']       = 1.0
           #namelistVals['R0']       = 1.0
           #namelistVals['RZ0']      = 0.0
            namelistVals['ZBOXMID']  = statedata['ZMID']    # 0.0
            namelistVals['ZBOXLEN']  = statedata['ZLEN']
            namelistVals['RBOXLFT']  = statedata['RLFT']
            namelistVals['RBOXLEN']  = statedata['RLEN']

            namelistVals['NCSCAL']   = 4
            namelistVals['PREDGE']   = statedata['pressure'][-1]

            importedVals['Jprl']     = statedata['jpar']
            importedVals['pressure'] = statedata['pressure']

        icntr  = 0
        niters = 25
        err_tol = 1.0e-6
        pre_err = 0.0

        t1 = time()
        print("CHEASE Start Time: ",datetime.fromtimestamp(t1))
        while True:
            cheasewrapper.init_chease_inputs(srcVals,namelistVals,importedVals)
            #--- run chease
            task_id = services.launch_task(1, self.cwd, chease_bin, logfile = 'xchease.log')
            retcode = services.wait_task(task_id)
            if (retcode != 0): raise Exception('Error executing chease')
            istatus = cheasewrapper.update_chease_outputs(suffix=icntr)

            if not init_run: break

            eqdskdata = cheasefiles.read_eqdsk(fpath=glob.glob('EQDSK_COCOS_??_POS.OUT')[0])
            ip_err = (importedVals['ITEXP'] - eqdskdata['CURNT'])/importedVals['ITEXP']

            if icntr >= niters:                    break
            elif abs(ip_err) <= err_tol:           break
            elif abs(pre_err - ip_err) <= 1.0e-6 : break
            else:
                pre_err = ip_err
                icntr += 1

                expeqdata = cheasefiles.read_expeq(fpath='EXPEQ.OUT.TOR')
                q         = cheasetools.interp(eqdskdata['rhopsi'],eqdskdata['q'],statedata['rho'])
                Jprl      = expeqdata['Jprl']*importedVals['B0EXP']/importedVals['R0EXP']/mu0
                pressure  = expeqdata['pressure']*importedVals['B0EXP']**2/mu0

                importedVals['ZMAX']     = expeqdata['zgeom'] * eqdskdata['RCTR']
                importedVals['Jprl']     = Jprl
                importedVals['pressure'] = pressure

                namelistVals['NCSCAL']   = 1 
                namelistVals['QSPEC']    = (1.0 - ip_err) * q[0]
                namelistVals['CSSPEC']   =  0.0

                if eq_recycle:
                    namelistVals['NOPT']     = -2 # 1
                    shutil.copyfile("NOUT","NIN")

        t2 = time()
        print("CHEASE End Time: ",datetime.fromtimestamp(t2))
        print("CHEASE Running Wall-Time = %3.3f" % (t2-t1))

        #--- update local geqdsk state
        cheaseqfname = glob.glob('EQDSK_COCOS_??_POS.OUT')[0]
        shutil.copyfile(cheaseqfname, self.eqdskfname)
        
        #--- load geqdsk to plasma state file
        do_update_state = int(getattr(self, "UPDATE_STATE", "1"))

        if do_update_state:
            ps = plasmastate('ips',1)
            ps.read(self.statefname)
            ps.load_geqdsk(self.eqdskfname)
           #ps.load_geqdsk(self.eqdskfname, keep_cur_profile=True, bdy_crat=1.e-6, kcur_option=1, rho_curbrk=0.9, EQ_Code_Info='chease')
            ps.store(self.statefname)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid,self.OUTPUT_FILES)

        #--- code exit
        print('cheasepy.step() done')

    def finalize(self, timeid=0):
        print('cheasepy.finalize() called')


