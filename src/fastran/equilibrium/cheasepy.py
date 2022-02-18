"""
 -----------------------------------------------------------------------
 chease component
 -----------------------------------------------------------------------
"""

import os
import sys
import shutil
import efit_io
import efit_eqdsk
import cheasefiles
import cheasewrapper

from numpy        import sqrt
from Namelist     import Namelist
from plasmastate  import plasmastate
from component import Component
#from ipsframework import Component

class cheasepy(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('cheasepy.init() started')
        services = self.services

        #--- initial force free equilibrium
        init_run = int(getattr(self, "INIT_RUN", 0))
        print('init_run = ', init_run)

        if init_run:
            print ('INIT CHEASEPY, force free')
            services.stage_state()
            cheasefiles.initial_equilibrium()
            services.update_state()

        #--- get shot and time
        self.TOKAMAK = services.get_config_param('TOKAMAK_ID')
        self.SHOT_ID = services.get_config_param('SHOT_NUMBER')
        self.TIME_ID = self.TIME_ID if hasattr(self, "TIME_ID") else services.get_config_param('TIME_ID')

        #--- use_instate
        self.useinstate       = int(getattr(self, "USE_INSTATE", "0"))
        if self.useinstate in [1,"YES"]: self.useinstate = 1
        else:                            self.useinstate = 0

        print('cheasepy.init() done')

    def step(self, timeid=0):
        print('cheasepy.step() started')
        #--- code entry
        services = self.services

        #--- GENE EXCUTABLE
        chease_bin = os.path.join(self.BIN_PATH, self.BIN)
        print("CHEASE Binary: %s" % chease_bin)

        #--- stage plasma state files
        services.stage_state()

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- get working directory
        self.cwd = services.get_working_dir()

        #--- get plasma state file names
        self.plasma_state_dir = services.get_config_param('PLASMA_STATE_WORK_DIR')
        if not self.useinstate:
           #self.eqdskfname       = services.get_config_param('CURRENT_EQDSK')
            self.statefname       = services.get_config_param('CURRENT_STATE')
           #self.fastranfname     = services.get_config_param('CURRENT_FASTRAN')
            self.boundaryfname    = services.get_config_param('CURRENT_BC')
           #self.path_to_eqdsk    = os.path.join(self.plasma_state_dir,self.eqdskfname)
            self.path_to_state    = os.path.join(self.plasma_state_dir,self.statefname)
           #self.path_to_fastran  = os.path.join(self.plasma_state_dir,self.fastranfname)
            self.path_to_boundary = os.path.join(self.plasma_state_dir,self.boundaryfname)
        else:
            self.instatefname     = services.get_config_param('CURRENT_INSTATE')
            self.path_to_instate  = os.path.join(self.plasma_state_dir,self.instatefname)
           #self.eqdskfname       = "g%s.%s" % (self.SHOT_ID,self.TIME_ID)

        #--- create inefit file
       #efit_io.io_input_from_state(self.statefname,self.boundaryfname)
       #inefit  = Namelist("inefit")
       #rho     = inefit["inefit"]["rho"]     # []
       #press   = inefit["inefit"]["press"]   #[Pa]
       #jpar    = inefit["inefit"]["jpar"]    #[A/m^2]

        #--- read inputs
        self.inputfpath = getattr(self,'INPUT_DIR',  '')
        self.inputfname = getattr(self,'INPUT_FILES','inchease').split()

        if 'inchease' in self.inputfname:
            self.incheasefname = 'inchease'
            self.path_to_inchease = os.path.join(self.inputfpath,self.incheasefname)
            if os.path.isfile(self.path_to_inchease):
               inchease = Namelist(self.path_to_inchease, case="lower")
            else:
               inchease = {}
        else:
            print("IOError: inchease file is missing. EXIT!"); sys.exit()

        if self.useinstate and 'instate' in self.inputfname:
            self.instatefname    = self.inputfname[self.inputfname.index('instate')]
            self.path_to_instate = os.path.join(self.inputfpath,self.instatefname)
        elif self.useinstate:
            print("IOError: instate file is missing. EXIT!"); sys.exit()

        #--- create input plasma profiles
        if not self.useinstate:
            statedata = cheasefiles.get_plasmastate(statefpath=self.path_to_instate,bcfpath=self.path_to_boundary)
           #self.iterdbfname = cheasefiles.update_iterdb_from_fastran(self.path_to_fastran,self.TOKAMAK,self.SHOT_ID,self.TIME_ID)
        elif   self.useinstate:
            statedata = cheasefiles.get_plasmastate(statefpath=self.path_to_instate)
           #statedata = cheasefiles.get_plasmastate(statefpath=self.path_to_instate,setParam={'mode':'non-kinetic'})

        namelistVals = {}
        if 'namelist' in inchease.keys():
            for ikey in inchease['namelist'].keys():
                namelistVals[ikey.upper()] = inchease['namelist'][ikey][0]

        #--- set default for namelist values
        self.NRBOX = int(float(getattr(self, "NRBOX", "129")))
        self.NZBOX = int(float(getattr(self, "NZBOX", "129")))
        self.RELAX = int(float(getattr(self, "RELAX", "0.0")))

        if   'NRBOX' not in namelistVals:         namelistVals['NRBOX'] = self.NRBOX
        elif namelistVals['NRBOX'] != self.NRBOX: namelistVals['NRBOX'] = self.NRBOX
        if   'NZBOX' not in namelistVals:         namelistVals['NZBOX'] = self.NZBOX
        elif namelistVals['NZBOX'] != self.NZBOX: namelistVals['NZBOX'] = self.NZBOX
        if   'RELAX' not in namelistVals:         namelistVals['RELAX'] = self.RELAX
        elif namelistVals['RELAX'] != self.RELAX: namelistVals['RELAX'] = self.RELAX

        srcVals = {}
        if 'sources'  in inchease.keys():
            for ikey in inchease['sources'].keys():
                srcVals[ikey] = inchease['sources'][ikey][0]

        #--- set default for source values
        self.rhomesh_src   = getattr(self,'rhomesh_src',    'eqdsk')
        self.current_src   = getattr(self,'current_src',    'eqdsk')
        self.pressure_src  = getattr(self,'pressure_src',   'eqdsk')
        self.boundary_src  = getattr(self,'boundary_src',   'eqdsk')
        self.boundary_type = getattr(self,'boundary_type',   'asis')
        self.eprofiles_src = getattr(self,'eprofiles_src', 'iterdb')
        self.iprofiles_src = getattr(self,'iprofiles_src', 'iterdb')

        if 'rhomesh_src'   not in srcVals: srcVals['rhomesh_src']   = self.rhomesh_src
        if 'current_src'   not in srcVals: srcVals['current_src']   = self.current_src
        if 'pressure_src'  not in srcVals: srcVals['pressure_src']  = self.pressure_src
        if 'boundary_src'  not in srcVals: srcVals['boundary_src']  = self.boundary_src
        if 'boundary_type' not in srcVals: srcVals['boundary_type'] = self.boundary_type
        if 'eprofiles_src' not in srcVals: srcVals['eprofiles_src'] = self.eprofiles_src
        if 'iprofiles_src' not in srcVals: srcVals['iprofiles_src'] = self.iprofiles_src

       #if not self.useinstate:
       #   srcVals['gfname']        = self.eqdskfname
       #   srcVals['inputpath']     = self.cwd
       #   srcVals['iterdbfname']   = self.iterdbfname

        importedVals = {}
        if 'imported'  in inchease.keys():
            for ikey in inchease['imported'].keys():
                importedVals[ikey] = inchease['imported'][ikey][0]

        if 'rhotor' in importedVals:
            namelistVals["NFUNRHO"]   = 1
            namelistVals["NRHOMESH"]  = 1
        else:
            namelistVals["NFUNRHO"]   = 0
            namelistVals["NRHOMESH"]  = 0

        importedVals['Te']       = statedata['Te']
        importedVals['ne']       = statedata['ne']
        importedVals['Ti']       = statedata['Ti']
        importedVals['ni']       = statedata['ni']
        importedVals['nz']       = statedata['nz']
        importedVals['Zeff']     = statedata['zeff']
        importedVals['Jprl']     = statedata['jpar']
        importedVals['ZMAX']     = statedata['zlim'][0]
        importedVals['ITEXP']    = statedata['ip']
        importedVals['R0EXP']    = statedata['RCTR']
        importedVals['B0EXP']    = statedata['BCTR']
        if 'rhotor' in importedVals:
           importedVals['rhotor']= statedata['rhotor']
        else:
           importedVals['rhopsi']= statedata['rhopsi']
        importedVals['rbound']   = statedata['rbound']
        importedVals['zbound']   = statedata['zbound']
        if 'ffprime' in statedata:
            importedVals['ffprime']  = statedata['ffprime']
        importedVals['pressure'] = statedata['pressure']

        cheasewrapper.init_chease_inputs(srcVals,namelistVals,importedVals)

        #--- run chease
        task_id = services.launch_task(1, self.cwd, chease_bin, logfile = 'xchease.log')
        retcode = services.wait_task(task_id)

        if (retcode != 0):
            raise Exception('Error executing chease')

        #--- update local geqdsk state
        shutil.copyfile('EQDSK_COCOS_02_POS.OUT', self.eqdskfname)
        
        #--- update local profile files
       #if not self.useinstate:
       #    self.path_to_iterdb = os.path.join(self.plasma_state_dir, "p%s.%s" % (self.SHOT_ID,self.TIME_ID))
       #    shutil.copyfile(self.iterdbfname, self.path_to_iterdb)

        #--- load geqdsk to plasma state file
        do_update_state = int(getattr(self, "UPDATE_STATE", "0"))

        if do_update_state:
            print('CHEASE: LOAD GEQDSK')
            ps = plasmastate('ips',1)
            ps.read(self.statefname)
            ps.load_geqdsk(self.eqdskfname)
            ps.store(self.statefname)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid,self.OUTPUT_FILES)

        #--- code exit
        print('cheasepy.step() done')

    def finalize(self, timeid=0):
        print('cheasepy.finalize() called')

   #def __call__(self,Component,**kwargs):
   #    services = self.services
   #    try:
   #        services.call('CHEASEQ', 'init', 0)
   #    except Exception:
   #        raise
   #   #self.init()
   #   #self.step()
   #   #self.finalize()
   #    return 1

#if __name__=='__main__':
#    cheasepy(Component,config)

