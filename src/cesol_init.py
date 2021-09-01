"""
 -----------------------------------------------------------------------
 CESOL init component
 -----------------------------------------------------------------------
"""

from component import Component

class cesol_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timestamp=0):
        print('>>> cesol_init.init() called')

    def step(self, timeStamp=0):
        services = self.services

        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')
        run_id = services.get_config_param('RUN_ID')

        print('tokamak_id =', tokamak_id)
        print('shot_number =',shot_number)
        print('run_id =', run_id)

        plasma_state_files = services.get_config_param('PLASMA_STATE_FILES')

        print('plasma state files =', plasma_state_files)

        self.INIT_PLASMA_STATE_FILES = getattr(self, 'INIT_PLASMA_STATE_FILES', '0')
        if self.INIT_PLASMA_STATE_FILES != 0:
            for plasma_state_file in plasma_state_files.split():
                print(plasma_state_file)
                open(plasma_state_file,"w").close()

            services.update_state()

    def checkpoint(self, timeStamp=0):
        print('>>> cesol_init.checkpoint() called')

        services = self.services
        services.stage_state()
        services.save_restart_files(timeStamp, self.RESTART_FILES)

    def finalize(self, timeStamp=0):
        print('>>> cesol_init.finalize() called')
