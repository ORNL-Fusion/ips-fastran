"""
 -----------------------------------------------------------------------
 fastran init component
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
from ipsframework import Component

class fastran_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timestamp=0.0):
        print('fastran_init.init() called')

    def step(self, timeStamp):
        #-- entry
        print('fastran_init.step() started')
        services = self.services

        #-- run identifiers
        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')

        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)

        #-- stage input files
        input_dir_id = getattr(self, "INPUT_DIR_ID", "")
        if input_dir_id:
            services.component_ref.config['INPUT_DIR'] = services.component_ref.config['INPUT_DIR']+"_%d"%int(float(input_dir_id))
        print('INPUT_DIR =', services.component_ref.config['INPUT_DIR'])

        services.stage_input_files(self.INPUT_FILES)

        #-- plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        #-- set default
        init_method = getattr(self, 'INIT_METHOD', 'instate') # not used
        instate_method = getattr(self, 'INSTATE_METHOD', '')
        f_instate = getattr(self, 'INSTATE', '')
        f_inps = getattr(self, 'INPS', '')
        f_ingeqdsk = getattr(self, 'INGEQDSK', '')
        f_eped = getattr(self, 'EPED', '')
        f_inpedestal = getattr(self, 'INPEDESTAL', '')

        #-- instate / dakota binding
        instate = Namelist(f_instate)

        var_list = [ var.upper() for var in instate["instate"].keys() ]
        for key in var_list:
            if hasattr(self, key):
                instate["instate"][key][0] = float(getattr(self, key))
                print(key, 'updated')

            if hasattr(self, "SCALE_"+key):
                scale = float(getattr(self, "SCALE_"+key))
                instate["instate"][key] = scale*np.array(instate["instate"][key])
                print(key,'scaled',scale)

        var_list = [ var.upper() for var in instate["flux"].keys() ]
        for key in var_list:
            if hasattr(self, "FLUX_"+key):
                instate["flux"][key][0] = float(getattr(self, "FLUX_"+key))
                print(key, 'updated')

        instate.write(f_instate)

        #-- expand instate
        #if f_inpedestal is not "":
        #    print('INPEDESTAL = ', f_inpedestal)
        #    update_instate_pedestal(f_instate, f_inpedestal)

        if instate_method == 'model':
            instate_model.instate_model(f_instate)

        #-- alloc plasma state file
        ps = plasmastate('ips', 1)

        instate = Namelist(f_instate)
        nicrf_src =  namelist_default(instate, "instate", "nicrf_src", [1])
        nlhrf_src =  namelist_default(instate, "instate", "nlhrf_src", [1])

        print('init from instate:', f_instate)

        ps.init(
            cur_state_file,
            global_label = ['ips'],
            runid = ['ips'],
            shot_number = [int(shot_number)],
            nspec_th = [instate["instate"]["n_ion"][0] + instate["instate"]["n_imp"][0]] ,
            nspec_beam = [instate["instate"]["n_beam"][0]],
            nspec_fusion = [instate["instate"]["n_fusion"][0]],
            nspec_rfmin = [instate["instate"]["n_min"][0]],
            nspec_gas = [instate["instate"]["n_ion"][0]],
            z_ion = instate["instate"]["z_ion"] + instate["instate"]["z_imp"],
            a_ion = instate["instate"]["a_ion"] + instate["instate"]["a_imp"],
            z_beam = instate["instate"]["z_beam"],
            a_beam = instate["instate"]["a_beam"],
            z_fusion = instate["instate"]["z_fusion"],
            a_fusion = instate["instate"]["a_fusion"],
            z_rfmin = instate["instate"]["z_min"],
            a_rfmin = instate["instate"]["a_min"],
            z_gas = instate["instate"]["z_ion"],
            a_gas = instate["instate"]["a_ion"],
            nicrf_src = nicrf_src,
            nlhrf_src = nlhrf_src,
            nrho = [101],
            time = [0.0],
            nlim = instate["instate"]["nlim"]
        )

        #-- input equlibrium
        if f_ingeqdsk:
            ps.load_geqdsk(f_ingeqdsk, keep_cur_profile=False, bdy_crat=0.01)
            shutil.copyfile(f_ingeqdsk, cur_eqdsk_file)

            geq = readg(f_ingeqdsk)
            r0  = geq["rzero" ]
            b0  = abs(geq["bcentr"])
            ip  = geq['cpasma']

            ps.load_j_parallel(instate["instate"]["rho"], instate["instate"]["j_tot"], "rho_eq", "curt", r0, b0, tot=True)
        else:
            r0 = instate["instate"]["r0"][0]
            b0 = instate["instate"]["b0"][0]

            rb = np.array(instate["instate"]["rbdry"])
            zb = np.array(instate["instate"]["zbdry"])

            R = 0.5*( np.max(rb) + np.min(rb) )
            Z = 0.5*( np.max(zb) + np.min(zb) )
            a = 0.5*( np.max(rb) - np.min(rb) )
            kappa = 0.5*( ( np.max(zb) - np.min(zb) )/ a )
            delta_u = ( R - rb[np.argmax(zb)] )/a
            delta_l = ( R - rb[np.argmin(zb)] )/a
            delta = 0.5*( delta_u + delta_l )

            ps.analytic_volume(b0, r0, a, kappa, delta)

            open(cur_eqdsk_file,"w").close()

        #-- load instate to plasma state
        print('instate to ps:', f_instate)
        instate_io.instate_to_ps(f_instate, ps)

        #-- write plasma state file
        ps.store(cur_state_file)

        #-- boundary condition state file
        instate = Namelist(f_instate)["inbc"]
        if not instate: instate = Namelist(f_instate)["instate"]

        inbc = Namelist()
        inbc["inbc"]["r0"] = instate["r0"]
        inbc["inbc"]["b0"] = instate["b0"]
        inbc["inbc"]["ip"] = instate["ip"]
        inbc["inbc"]["nbdry"] = instate["nbdry"]
        inbc["inbc"]["rbdry"] = instate["rbdry"]
        inbc["inbc"]["zbdry"] = instate["zbdry"]
        inbc["inbc"]["nlim"] = instate["nlim"]
        inbc["inbc"]["rlim"] = instate["rlim"]
        inbc["inbc"]["zlim"] = instate["zlim"]
        inbc.write(cur_bc_file)

        #-- instate, fastran nc file
        if f_instate: shutil.copyfile(f_instate, cur_instate_file)

        #-- eped
        if f_eped:
            cesol_io.eped_to_plasmastate(cur_state_file, f_eped, cur_instate_file)

        #-- update plasma state
        services.update_state()

        #-- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0.0):
        print('fastran_init.finalize() called')
