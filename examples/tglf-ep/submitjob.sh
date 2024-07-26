#!/bin/bash -l
module load python
source activate /global/common/software/atom/perlmutter/cesol/conda/latest
export PYTHONPATH=/pscratch/sd/p/parkjm/tglf-ep/ips-fastran/src:$PYTHONPATH

# export TGLFEP_BIN_DIR='/global/common/software/atom/gacode_add_TGLFEP/TGLF-EP'
# export TGLFEP_BIN_NAME='TGLFEP_driver'
export TGLFEP_BIN_DIR='/pscratch/sd/p/parkjm/tglf-ep/ips-fastran/examples/tglf-ep'
export TGLFEP_BIN_NAME='TGLFEP_driver'

export SHOT_NUMBER=000001
export TIME_ID=00001

ips.py --simulation=fastran_scenario.config --platform=perlmutter_cpu_node.conf --log=ips.log 1> ips.out 2> ips.err &

conda deactivate
