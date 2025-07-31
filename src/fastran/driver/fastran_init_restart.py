"""
 -----------------------------------------------------------------------
 fastran init component, restart from plsmas state and geqdsk
 -----------------------------------------------------------------------
"""

import shutil
from Namelist import Namelist
from fastran.equilibrium.efit_eqdsk import readg
from ipsframework import Component


class fastran_init_restart (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print ('fastran_init_restart.init() called')

    def step(self, timeid=0):
        print ('fastran_init_restart.step() start')

        # -- run identifiers
        tokamak_id = self.services.get_config_param('TOKAMAK_ID')
        shot_number = self.services.get_config_param('SHOT_NUMBER')

        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        # -- plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = self.services.get_config_param('CURRENT_BC')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        try:
            cur_fastran_file = self.services.get_config_param('CURRENT_FASTRAN')
        except:
            cur_fastran_file = ''

        # -- set default
        f_inps = getattr(self, 'INPS', '')
        f_ingeqdsk = getattr(self, 'INGEQDSK', '')
        f_instate = getattr(self, 'INSTATE', '')

        if f_inps: shutil.copyfile(f_inps, cur_state_file)
        if f_ingeqdsk: shutil.copyfile(f_ingeqdsk, cur_eqdsk_file)
        if f_instate: shutil.copyfile(f_instate, cur_instate_file)
        if cur_fastran_file: open(cur_fastran_file, "w").close()

        # -- boundary condition state file
        geq = readg(f_ingeqdsk)
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']

        inbc = Namelist()
        inbc["inbc"]["r0"] = [r0]
        inbc["inbc"]["b0"] = [b0]
        inbc["inbc"]["ip"] = [ip*1.0e-6]
        inbc["inbc"]["nbdry"] = [geq["nbdry"]]
        inbc["inbc"]["rbdry"] = geq["rbdry"]
        inbc["inbc"]["zbdry"] = geq["zbdry"]
        inbc["inbc"]["nlim"] = [geq["nlim"]]
        inbc["inbc"]["rlim"] = geq["rlim"]
        inbc["inbc"]["zlim"] = geq["zlim"]
        inbc.write(cur_bc_file)

        # -- update plasma state
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

        # -- exit
        print('fastran_init_restart.step() done')

    def finalize(self, timeid=0):
        print('fastran_init_restart.finalize() called')
