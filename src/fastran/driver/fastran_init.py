"""
 -----------------------------------------------------------------------
 fastran init component
 -----------------------------------------------------------------------
"""

import shutil
import numpy as np
import netCDF4
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.equilibrium.efit_eqdsk import readg
from fastran.instate import instate_model
from fastran.driver import cesol_io
from fastran.util.fastranutil import namelist_default
from fastran.util import dakota_io
from fastran.state.instate import Instate


class fastran_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('fastran_init.init() called')
        state_files = self.services.get_config_param('STATE_FILES')
        print(state_files)
        for fname in state_files.split():
            print(f'create null {fname}')
            open(fname, 'w').close()

    def step(self, timeid=0):
        print('fastran_init.step() started')

        # -- run identifiers
        tokamak_id = self.services.get_config_param('TOKAMAK_ID')
        shot_number = self.services.get_config_param('SHOT_NUMBER')

        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)

        # -- stage input files
        input_dir_id = getattr(self, 'INPUT_DIR_ID', '')
        if input_dir_id:
            self.services.component_ref.config['INPUT_DIR'] = self.services.component_ref.config['INPUT_DIR'] + \
                '_%d' % int(float(input_dir_id))
        print('INPUT_DIR =', self.services.component_ref.config['INPUT_DIR'])

        self.services.stage_input_files(self.INPUT_FILES)

        # -- plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        # cur_bc_file = self.services.get_config_param('CURRENT_BC')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        # -- set default
        init_method = getattr(self, 'INIT_METHOD', 'instate')  # not used
        instate_method = getattr(self, 'INSTATE_METHOD', '')
        f_instate = getattr(self, 'INSTATE', '')
        f_inps = getattr(self, 'INPS', '')
        f_ingeqdsk = getattr(self, 'INGEQDSK', '')
        f_eped = getattr(self, 'EPED', '')
        f_inpedestal = getattr(self, 'INPEDESTAL', '')

        # -- instate / dakota binding
        instate = Instate(f_instate)

        print("start dakota update")
        dakota_io.update_namelist(self, instate.data, section="instate")
        dakota_io.update_namelist(self, instate.data)

        # -- expand the instate profiles
        if instate_method == 'model':
            print("instate_method: model")
            instate.set_shape()
            instate.model_profile()
        elif instate_method == 'trace':
            print("instate_method: trace")
            tmin = float(self.services.sim_conf['ITERATION_LOOP']['TMIN'])
            tmax = float(self.services.sim_conf['ITERATION_LOOP']['TMAX'])
            dt = float(self.services.sim_conf['ITERATION_LOOP']['DT'])
            f_timetrace = getattr(self, 'TIMETRACE', 'timetrace.nc')
            interpolation_method = getattr(self, 'INTERPOLATION_METHOD', 'linear')
            timenow = tmin
            instate.set_grid()
            instate.zeros()
            instate.from_timetrace(f_timetrace, timenow, interpolation_method=interpolation_method, init=True)
            #instate.particle_balance()
        else:
            instate.zeros()
            instate.write('instate_test')

        instate.particle_balance()
        instate.zeros()
        instate.scale()

        # -- alloc plasma state file
        ps = plasmastate('ips', 1)

        nicrf_src = namelist_default(instate.data, "instate", "nicrf_src", [1])
        nlhrf_src = namelist_default(instate.data, "instate", "nlhrf_src", [1])

        print('init from instate:', f_instate)
        ps.init(
            cur_state_file,
            global_label=['ips'],
            runid=['ips'],
            shot_number=[int(shot_number)],
            nspec_th=[instate["n_ion"][0] + instate["n_imp"][0]],
            nspec_beam=[instate["n_beam"][0]],
            nspec_fusion=[instate["n_fusion"][0]],
            nspec_rfmin=[instate["n_min"][0]],
            nspec_gas=[instate["n_ion"][0]],
            z_ion=np.append(instate["z_ion"], instate["z_imp"]),
            a_ion=np.append(instate["a_ion"], instate["a_imp"]),
            z_beam=instate["z_beam"],
            a_beam=instate["a_beam"],
            z_fusion=instate["z_fusion"],
            a_fusion=instate["a_fusion"],
            z_rfmin=instate["z_min"],
            a_rfmin=instate["a_min"],
            z_gas=instate["z_ion"],
            a_gas=instate["a_ion"],
            nicrf_src=nicrf_src,
            nlhrf_src=nlhrf_src,
            nrho=[101],
            time=[0.0],
            nlim=instate["nlim"]
        )

        # -- input equlibrium
        if f_ingeqdsk:
            ps.load_geqdsk(f_ingeqdsk, keep_cur_profile=False, bdy_crat=0.01)
            shutil.copyfile(f_ingeqdsk, cur_eqdsk_file)

            geq = readg(f_ingeqdsk)
            r0 = geq["rzero"]
            b0 = abs(geq["bcentr"])
            ip = geq['cpasma']

            ps.load_j_parallel(instate["rho"], instate["j_tot"], "rho_eq", "curt", r0, b0, tot=True)
        else:
            r0 = instate["r0"][0]
            b0 = instate["b0"][0]

            rb = np.array(instate["rbdry"])
            zb = np.array(instate["zbdry"])

            R = 0.5*(np.max(rb) + np.min(rb))
            Z = 0.5*(np.max(zb) + np.min(zb))
            a = 0.5*(np.max(rb) - np.min(rb))
            kappa = 0.5*((np.max(zb) - np.min(zb)) / a)
            delta_u = (R - rb[np.argmax(zb)])/a
            delta_l = (R - rb[np.argmin(zb)])/a
            delta = 0.5*(delta_u + delta_l)

            ps.analytic_volume(b0, r0, a, kappa, delta)

            open(cur_eqdsk_file, "w").close()
            print('end')

        # -- load instate to plasma state
        # instate.model_beam = 1
        instate.to_ps(ps, init=True)

        # -- write plasma state file
        ps.store(cur_state_file)

        # -- recycling neutral
        neutral_energy = namelist_default(instate.data, "instate", "neutral_energy", [15.0])[0]
        neutral_recycling = namelist_default(instate.data, "instate", "neutral_recycling", [1.0e20])[0]

        ps = netCDF4.Dataset(cur_state_file, "r+", format='NETCDF4')
        ps.variables["sc0_to_sgas"][:] = [1]
        ps.variables["is_recycling"][:] = [1]
        ps.variables["sc0"][:] = [neutral_recycling]
        ps.variables["e0_av"][:] = [neutral_energy*1.0e-3]
        ps.close()

        # -- write instate file
        instate.write(f_instate)

        # -- instate, fastran nc file
        if f_instate:
            shutil.copyfile(f_instate, cur_instate_file)

        #-- eped
        if f_eped:
            cesol_io.eped_to_plasmastate(cur_state_file, f_eped, cur_instate_file)

        # -- update plasma state
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid):
        print('fastran_init.finalize() called')
