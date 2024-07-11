"""
 -----------------------------------------------------------------------
 fastran timetrace init component
 -----------------------------------------------------------------------
"""

import shutil
import numpy as np
from Namelist import Namelist
from ipsframework import Component
from fastran.plasmastate.plasmastate import plasmastate
from fastran.util.fastranutil import namelist_default

from fastran.state.instate import Instate
from fastran.state.timetrace import Timetrace


class timetrace_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> timetrace_init.init() called')

    def step(self, timeid=0):
        # -- entry
        print('>>> timetrace_init.step() started')

        # -- run identifiers
        tokamak_id = self.services.get_config_param('TOKAMAK_ID')
        shot_number = self.services.get_config_param('SHOT_NUMBER')
        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)

        self.services.stage_input_files(self.INPUT_FILES)

        # -- plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        next_state_file = self.services.get_config_param('NEXT_STATE')
        next_eqdsk_file = self.services.get_config_param('NEXT_EQDSK')
        next_instate_file = self.services.get_config_param('NEXT_INSTATE')

        # -- instate / dakota binding
        f_instate = getattr(self, 'INSTATE', '')
        if f_instate and f_instate is not cur_instate_file:
            shutil.copyfile(f_instate, cur_instate_file)

        # -- alloc plasma state file
        instate = Instate(cur_instate_file)
        instate.set_grid()
#       instate.model_profile()
        instate.zeros()
#       instate.particle_balance()

        print('init from instate:', f_instate)
        ps = plasmastate('ips', 1)
        nicrf_src = namelist_default(instate.data, "instate", "nicrf_src", [1])
        nlhrf_src = namelist_default(instate.data, "instate", "nlhrf_src", [1])
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

        # -- from timetrace
        f_timetrace = self.TIMETRACE
        timenow = float(self.services.sim_conf['ITERATION_LOOP']['TMIN'])
        instate.from_timetrace(
            f_timetrace,
            timenow,
            interpolation_method='linear')
        instate.particle_balance()

        # -- analytic volume
        rb = np.array(instate['rbdry'])
        zb = np.array(instate['zbdry'])

        R = 0.5 * (np.max(rb) + np.min(rb))
        Z = 0.5 * (np.max(zb) + np.min(zb))
        a = 0.5 * (np.max(rb) - np.min(rb))
        kappa = 0.5 * ((np.max(zb) - np.min(zb)) / a)
        delta_u = (R - rb[np.argmax(zb)]) / a
        delta_l = (R - rb[np.argmin(zb)]) / a
        delta = 0.5 * (delta_u + delta_l)

        b0 = instate['b0'][0]
        r0 = instate['r0'][0]

        ps.analytic_volume(b0, r0, a, kappa, delta)

        # -- instate to ps
        instate.to_ps(ps)

        # -- write plasma state file
        instate.write(cur_instate_file)
        ps.store(cur_state_file)

        # -- touch plasma state files
        for fname in [
                cur_eqdsk_file,
                next_state_file,
                next_eqdsk_file,
                next_instate_file]:
            open(fname, 'w').close()

        # -- update plasma state
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('>>> timetrace_init.finalize() called')
