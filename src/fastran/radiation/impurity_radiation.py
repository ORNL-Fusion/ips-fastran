
'''
 -----------------------------------------------------------------------
 radiation component
 -----------------------------------------------------------------------
'''

import numpy as np
from scipy import optimize
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate
from fastran.util import dakota_io
from fastran.radiation.radiation_io import radiation_io


class impurity_radiation(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('radiation.init() started')

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        # -- RADAS data directory
        self.radas_dir = self.ADAS_DATA_DIR

        # -- read ADAS data
        instate = Instate(cur_instate_file)
        a_imp = instate['a_imp']
        z_imp = instate['z_imp']
        spec_map = {2:'helium', 4:'beryllium', 6:'carbon', 7:'nitrogen', 10:'neon', 18:'argon', 36:'krypton', 74:'tungsten'}
        self.Lz = {}
        print(z_imp)
        for z in z_imp:
            spec = spec_map[z]
            print(f'z_imp = {z}, spec = {spec}')
            self.Lz[z] = radiation_io(self.radas_dir)
            self.Lz[z].read(spec)

    def step(self, timeid=0):
        print('radiation.step() started')

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- impurity radiation loss
        self.update_state(cur_state_file, cur_instate_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('radiation.finalize() called')

    def update_state(self, cur_state_file, cur_instate_file):
        instate = Instate(cur_instate_file)
        a_imp = instate['a_imp']
        z_imp = instate['z_imp']
        rho = instate['rho']
        for k, z in enumerate(z_imp):
            ne = instate[f'ne'] * 1.e19
            nz = instate[f'density_imp_{k}'] * 1.e19
            te = instate['te'] * 1.e3
            p_rad = self.Lz[z](te, ne) * ne * nz
            print(k, p_rad)
        instate['p_rad'] = -p_rad * 1.e-6
        instate.write(cur_instate_file)

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)
        ps.load_vol_profile(rho, p_rad, 'rho_rad', 'prad')
        ps.store(cur_state_file)