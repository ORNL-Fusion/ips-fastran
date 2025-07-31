"""
 -----------------------------------------------------------------------
 CESOL init component
 -----------------------------------------------------------------------
"""

from ipsframework import Component

class cesol_init (Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timestamp=0):
        print('>>> cesol_init.init() called')

    def step(self, timeid=0):
        print('>>> cesol_init.step() called')
        tokamak_id = self.services.get_config_param('TOKAMAK_ID')
        shot_number = self.services.get_config_param('SHOT_NUMBER')
        run_id = self.services.get_config_param('RUN_ID')

        print('tokamak_id =', tokamak_id)
        print('shot_number =', shot_number)
        print('run_id =', run_id)

        plasma_state_files = self.services.get_config_param('STATE_FILES')

        print('plasma state files =', plasma_state_files)

        self.INIT_PLASMA_STATE_FILES = getattr(self, 'INIT_PLASMA_STATE_FILES', '0')
        if self.INIT_PLASMA_STATE_FILES != 0:
            for plasma_state_file in plasma_state_files.split():
                print('touch '+plasma_state_file)
                open(plasma_state_file,"w").close()

            self.services.update_state()

    def checkpoint(self, timeid=0):
        print('>>> cesol_init.checkpoint() called')

        self.services.stage_state()
        self.services.save_restart_files(timeStamp, self.RESTART_FILES)

    def finalize(self, timeid=0):
        print('>>> cesol_init.finalize() called')
