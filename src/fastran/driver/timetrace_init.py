"""
 -----------------------------------------------------------------------
 fastran timetrace init component
 -----------------------------------------------------------------------
"""

import shutil
import numpy as np
from Namelist import Namelist
from fastran.plasmastate.plasmastate import plasmastate
from fastran.equilibrium.efit_eqdsk import readg
from fastran.instate import instate_model
from fastran.instate import instate_io
from fastran.driver import cesol_io
#from fastran.stability.pedestal_io import update_instate_pedestal
from fastran.util.fastranutil import namelist_default

from fastran.state.instate import Instate
from fastran.state.timetrace import Timetrace

from ipsframework import Component

class timetrace_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> timetrace_init.init() called')

    def step(self, timeid=0):
        #-- entry
        print('>>> timetrace_init.step() started')

        #-- run identifiers
        tokamak_id = self.services.get_config_param('TOKAMAK_ID')
        shot_number = self.services.get_config_param('SHOT_NUMBER')
        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)

        tmin = int(self.services.sim_conf["TIME_LOOP"]["TMIN"])
        tmax = int(self.services.sim_conf["TIME_LOOP"]["TMAX"])
        print('tmin =', tmin)
        print('tmax =', tmax)

        self.services.stage_input_files(self.INPUT_FILES)

        #-- plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_coil_file = self.services.get_config_param('CURRENT_COIL')

        #-- instate / dakota binding
        f_instate = getattr(self, 'INSTATE', '')
        if f_instate: shutil.copyfile(f_instate, cur_instate_file)

        #-- alloc plasma state file
        instate = Instate(cur_instate_file)
        instate.model_profile()
        instate.zeros()
#       instate.particle_balance()

        ps = plasmastate('ips', 1)
        ps.init(
            ps,
            global_label = ['ips'],
            runid = ['ips'],
            shot_number = [123456],
            nspec_th = [instate["n_ion"][0] + instate["n_imp"][0]] ,
            nspec_beam = [instate["n_beam"][0]],
            nspec_fusion = [instate["n_fusion"][0]],
            nspec_rfmin = [instate["n_min"][0]],
            nspec_gas = [instate["n_ion"][0]],
            z_ion = instate["z_ion"] + instate["z_imp"],
            a_ion = instate["a_ion"] + instate["a_imp"],
            z_beam = instate["z_beam"],
            a_beam = instate["a_beam"],
            z_fusion = instate["z_fusion"],
            a_fusion = instate["a_fusion"],
            z_rfmin = instate["z_min"],
            a_rfmin = instate["a_min"],
            z_gas = instate["z_ion"],
            a_gas = instate["a_ion"],
            nicrf_src = instate["nicrf_src"],
            nlhrf_src = instate["nlhrf_src"],
            nrho = [101],
            time = [0.0],
            nlim = instate["nlim"]
        )

        #-- from timetrace 
        f_timetrace = self.TIMETRACE
        timetrace = Timetrace(f_timetrace, start_time=tmin, end_time=tmax)

        ip = timetrace.get("ip", timeid)
        b0 = timetrace.get("bt", timeid)
        r0 = timetrace.get("r0", timeid)
        print('ip = {}, bt={}, r0={}'.format(ip, b0, r0))

        instate["ip"] = [ip*1.e-6]
        instate["b0"] = [b0]
        instate["r0"] = [r0]
        instate["ne"] = timetrace.slice("ne", timeid)
        instate["te"] = timetrace.slice("te", timeid)
        instate["ti"] = timetrace.slice("ti", timeid)
        instate["omega"] = timetrace.slice("omega", timeid)
        instate["zeff"] = 101*[1.6]

        instate.particle_balance()

        instate["rbdry"] = timetrace.get("rbdry", timeid)
        instate["zbdry"] = timetrace.get("zbdry", timeid)
        instate["nbdry"] = [timetrace.get("nbdry", timeid)]

        instate["p_eq"] = timetrace.get("p_eq", timeid)
        instate["j_tot"] = timetrace.get("jpar", timeid)*1.0e-6

        rb = np.array(instate["rbdry"])
        zb = np.array(instate["zbdry"])
        
        R = 0.5*( np.max(rb) + np.min(rb) )
        Z = 0.5*( np.max(zb) + np.min(zb) )
        a = 0.5*( np.max(rb) - np.min(rb) )
        kappa = 0.5*( ( np.max(zb) - np.min(zb) )/ a )
        delta_u = ( R - rb[np.argmax(zb)] )/a
        delta_l = ( R - rb[np.argmin(zb)] )/a
        delta = 0.5*( delta_u + delta_l )
        
        ps.analytic_volume(b0, r0, a, kappa, delta)

        instate.to_ps(ps)

        #-- write plasma state file
        instate.write(cur_instate_file)
        ps.store(cur_state_file)

        #-- instate, fastran nc file
        #if f_instate: shutil.copyfile(f_instate, cur_instate_file)

        #-- touch plasma state files
        for fname in [cur_eqdsk_file]:
            open(fname, "w").close()

        #-- touch plasma state files
        for fname in [cur_coil_file]:
            open(fname, "w").close()

        #-- update plasma state
        self.services.update_state()

        #-- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('>>> timetrace_init.finalize() called')
