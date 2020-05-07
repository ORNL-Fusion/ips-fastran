"""
 -----------------------------------------------------------------------
 dummy init component
 -----------------------------------------------------------------------
"""
from component import Component

class dummy_init (Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid):
        print('dummy_init.init() called')

    def step(self, timeid):
        print('dummy_init.step() started')

        services = self.services

        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')
        run_id = services.get_config_param('RUN_ID')

        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)
        print('run_id =', run_id)

        plasma_state_files = services.get_config_param('PLASMA_STATE_FILES')

        print('plasma state files =', plasma_state_files)

        self.INIT_PLASMA_STATE_FILES = getattr(self,'INIT_PLASMA_STATE_FILES', '0')
        if self.INIT_PLASMA_STATE_FILES != 0:
            for plasma_state_file in plasma_state_files.split():
                print(plasma_state_file)
                open(plasma_state_file, "w").close()

            services.update_plasma_state()

    def checkpoint(self, timeid):
        print('dummy_init.checkpoint() called')

        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timeid, self.RESTART_FILES)

    def finalize(self, timeid):
        print('dummy_init.finalize() called')
