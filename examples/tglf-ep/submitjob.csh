#!/bin/bash -l
module load python
source activate /global/common/software/atom/perlmutter/cesol/conda/latest
setenv PYTHONPATH /pscratch/sd/p/parkjm/tglf-ep/ips-fastran/src:$PYTHONPATH

# export TGLFEP_BIN_DIR='/global/common/software/atom/gacode_add_TGLFEP/TGLF-EP'
# export TGLFEP_BIN_NAME='TGLFEP_driver'
setenv TGLFEP_BIN_DIR '/pscratch/sd/p/parkjm/tglf-ep/ips-fastran/examples/tglf-ep'
setenv TGLFEP_BIN_NAME 'TGLFEP_driver'

setenv SHOT_NUMBER 000001
setenv TIME_ID 00001

ips.py --simulation=fastran_scenario.config --platform=perlmutter_cpu_node.conf --log=ips.log 1> ips.out 2> ips.err &

conda deactivate
