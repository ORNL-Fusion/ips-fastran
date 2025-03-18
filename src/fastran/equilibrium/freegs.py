"""
 -----------------------------------------------------------------------
 freegs component
 -----------------------------------------------------------------------
"""
import os
import shutil
import time as timer
import numpy as np
from Namelist import Namelist

from ipsframework import Component
from fastran.equilibrium import efit_io
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver.inmetric_io import ps_to_inmetric
from fastran.equilibrium.efit_eqdsk import readg
from fastran.util.fastranutil import freeze
from fastran.equilibrium.freegs_io import freegs_io
from fastran.state.instate import Instate


class freegs(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def _get_shot_time(self):
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(float(self.TIME_ID)) if hasattr(self, 'TIME_ID') else int(self.services.get_config_param('TIME_ID'))
        return ishot, itime

    def init(self, timeid=0):
        print('>>> freegs.init() started')
        self.ps_backend = getattr(self, 'PS_BACKEND', 'PS')

        # -- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid
        print('freegs time = {}'.format(itime))

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        #--- stage plasma state files
        self.services.stage_state()

        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')

        instate = Instate(cur_instate_file)
        if self.ps_backend == 'PS':
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            instate.from_ps(ps)
            instate.write(cur_instate_file)

        # -- init freegs
        rmin = float(getattr(self, 'RMIN', '0'))
        rmax = float(getattr(self, 'RMAX', '0'))
        zmax = float(getattr(self, 'ZMAX', '0'))
        zmin = float(getattr(self, 'ZMIN', '0'))
        domain = {'rmin':rmin, 'rmax':rmax, 'zmin':zmin, 'zmax':zmax}
        domain_valid = rmin > 0 and rmax > rmin and zmax - zmin > 0 and zmax > zmin
        print(f'domain = {domain}')
        if not domain_valid:
           raise Exception('freegs computational domain error')

        topology = getattr(self, 'TOPOLOGY', '')
        print(f'topology = {topology}')
        if topology not in ['DN', 'LSN', 'USN']:
           raise Exception('freegs topology error')

        nx = int(getattr(self, 'NX', '129'))
        ny = int(getattr(self, 'NY', '129'))

        incoil = getattr(self, 'INCOIL', 'coils.json')
        print(f'incoil = {incoil}')

        model_pressure = getattr(self, 'PRESSURE', 'kinetic')
        print(f'model_pressure = {model_pressure}')

        self.freegs_driver = freegs_io(topology=topology, model_pressure=model_pressure)
        self.freegs_driver.from_instate(cur_instate_file)
        self.freegs_driver.define_tokamak(f_coil_data=incoil, domain=domain, nx=nx, ny=ny)

        # -- initial equilibrium
        init_run = int(getattr(self, 'INIT_RUN', 0))
        print(f'init_run = {init_run}')
        if init_run:
            print('initial approximate equilibrium')
            self.freegs_driver.set_profiles_simple()
            self.freegs_driver.solve(p0=[0.1, 0.5, 1.0], verbose=True)
            self.freegs_driver.write_geqdsk(cur_eqdsk_file)
            if self.ps_backend == 'PS':
                ps = plasmastate('ips',1)
                ps.read(cur_state_file)
                ps.load_geqdsk(cur_eqdsk_file)
                ps.store(cur_state_file)
            elif self.ps_backend == 'INSTATE':
                instate.inmetric(cur_eqdsk_file)
            self.services.update_state()

    def step(self, timeid=0):
        "nrho = 257 (instead of 101) for PS backend"
        print('>>> freegs.step() started')

        # -- freeze/resume
        if freeze(self, timeid, 'freegs'): return None

        # -- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0: itime = timeid
        print('freegs time = {}'.format(itime))

        # -- stage plasma state files
        self.services.stage_state()

        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')

        instate = Instate(cur_instate_file)
        if self.ps_backend == 'PS':
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            instate.from_ps(ps)
            instate.write(cur_instate_file)

        # -- call freegs
        niter = int(getattr(self, 'NITER', '5'))
        print('niter = ', niter)

        self.freegs_driver.from_instate(cur_instate_file)
        for k in range(niter):
            relax = 0.5 if k > 0 else 0
            print(f'>>> freegs iteration {k}')
            self.freegs_driver.set_profiles(cur_eqdsk_file, relax=relax)
            self.freegs_driver.solve(method=1)
            self.freegs_driver.write_geqdsk(cur_eqdsk_file)

        # -- load geqdsk to state
        if self.ps_backend == 'PS':
            ps = plasmastate('ips',1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)
        elif self.ps_backend == 'INSTATE':
            instate.inmetric(cur_eqdsk_file)
            instate.write(cur_instate_file)
        self.services.update_state()

        # -- update state files
        self.services.update_state()

        #--- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid):
        print ('>>> freegs.finalize() called')
