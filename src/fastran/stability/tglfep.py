"""
 -----------------------------------------------------------------------
 tglfep component
 -----------------------------------------------------------------------
"""

import os
import shutil
from ipsframework import Component
from fastran.stability import tglfep_io


class tglfep(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('tglfep.init() called')

    def step(self, timeid=0):
        print('tglfep.step() started')

        # -- excutable
        tglfep_bin = os.path.join(self.BIN_PATH, self.BIN_TGLFEP)
        alpha_bin = os.path.join(self.BIN_PATH, self.BIN_ALPHA)
        print(tglfep_bin)
        print(alpha_bin)

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')

        # -- stage input files
        self.services.stage_input_files(self.INPUT_FILES)

        tglfep_io.write_inputfiles(cur_state_file, cur_eqdsk_file)

        # -- run tglfep
        print('run tglfep')

        cwd = self.services.get_working_dir()
        task_id = self.services.launch_task(self.NPROC, cwd, tglfep_bin, logfile='tglfep.log')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            print(retcode)
            raise Exception('Error executing: tglfep')
        task_id = self.services.launch_task(self.NPROC, cwd, alpha_bin, logfile='alpha.log')
        retcode = self.services.wait_task(task_id)
        if (retcode != 0):
            print(retcode)
            raise Exception('Error executing: Alpha')

        # -- get tglfep output
#        tglfep_io.update_state()
        from Namelist import Namelist

        print('tglfep update_state')
    
        # read tglfep output

        with open('alpha_dpdr_crit.input', 'r') as f:
            cg_str = f.readlines()
            f.close()

        nr = len(cg_str) - 1
        cg_list = [0.] * nr
        for ir in range(nr):
            cg_list[ir] = cg_str[ir+1]

        with open('density_alpha.out', 'r') as f:
            den_str = f.readlines()
            f.close()

        nr = len(den_str) - 1
        den_list = [0.] * nr
        for ir in range(nr):
            den_list[ir] = den_str[ir+1]

        with open('pe_fus.out', 'r') as f:
            pe_str = f.readlines()
            f.close()

        nr = len(pe_str) - 1
        pe_list = [0.] * nr
        for ir in range(nr):
            pe_list[ir] = pe_str[ir+1]

        with open('pi_fus.out', 'r') as f:
            pi_str = f.readlines()
            f.close()

        nr = len(pi_str) - 1
        pi_list = [0.] * nr
        for ir in range(nr):
            pi_list[ir] = pi_str[ir+1]

#        with open('alpha_flow.out', 'r') as f:
#            flow_str = f.readlines()
#            f.close
#
#        nr = len(flow_str) - 1
#        flow_list = [0.] * nr
#        for ir in range(nr):
#            flow_list[ir] = flow_str[ir+1]

        # update state file

        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        instate = Namelist(cur_instate_file)
        instate['EP']['alpha_critical_gradient'] = cg_list
        instate['EP']['critical_gradient_units'] = '10 kPa/m'
#        instate['EP']['alpha_flow'] = flow_list
#        instate['EP']['alpha_flow_units'] = '10^19/s'
        instate['pe_fus'] = pe_list
        instate['pi_fus'] = pi_list
        instate['density_alpha'] = den_list

        instate.write(cur_instate_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

    def finalize(self, timeid=0):
        print('tglfep.finalize() called')
