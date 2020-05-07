# nested workflow wrapper for fastran workflow

import os
import glob
import shutil

from  component import Component
import cesol_io

class call_fastran_driver(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('>>> entry: call_fastran_driver init')

        override = {}
        override['PWD'] = self.services.get_config_param('PWD')
        override['TOKAMAK_ID'] = self.services.get_config_param('TOKAMAK_ID')
        override['SHOT_NUMBER'] = self.services.get_config_param('SHOT_NUMBER')
        override['TIME_ID'] = self.services.get_config_param('TIME_ID')
        override['PLASMA_STATE_WORK_DIR'] = self.services.get_config_param('PLASMA_STATE_WORK_DIR')
        print(override)

        input_dir = 'input'
        os.mkdir(os.path.join(os.getcwd(),input_dir))

        input_dir_sim = self.services.get_config_param('INPUT_DIR_SIM')
        print('INPUT_DIR_SIM =', input_dir_sim)
        input_files = glob.glob(os.path.join(input_dir_sim,"*"))
        for f in input_files:
            print(f)
            shutil.copyfile(f,os.path.join(os.path.join(os.getcwd(), input_dir), os.path.basename(f)))

        (sim_name, init, driver) = \
            self.services.create_sub_workflow('fastran_driver', self.SUB_WORKFLOW, override, input_dir)

        self.init = init
        self.driver = driver
        print('sub workflow name:', sim_name)

        self.services.call(self.init, 'init', timeid)
        self.services.call(self.init, 'step', timeid)
        self.services.call(self.init, 'finalize', timeid)

        self.services.call(self.driver,'init', timeid)

    def step(self, timeid=0):
        print('>>> entry: call_fastra_driver step')

        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_fastran_file = self.services.get_config_param('CURRENT_FASTRAN')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eped_state = self.services.get_config_param('EPED_STATE')

        self.services.call(self.driver,'step', timeid)

        self.services.stage_plasma_state()
        cesol_io.fastran_to_instate(cur_fastran_file, cur_instate_file)
        self.services.update_plasma_state()

        self.services.stage_subflow_output_files()

    def finalize(self, timeid=0):
        print('>>> entry: call_fastra_driver finalize')

        self.services.call(self.driver, 'finalize', '0.0')
