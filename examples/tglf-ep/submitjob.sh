#!/bin/bash -l
source /usr/share/lmod/lmod/init/bash

module load python
source activate /global/common/software/atom/perlmutter/cesol/conda/latest
export PYTHONPATH=/global/homes/b/bassem/ips-fastran/src:$PYTHONPATH
echo $PYTHONPATH

# export TGLFEP_BIN_DIR='/global/common/software/atom/gacode_add_TGLFEP/TGLF-EP'
# export TGLFEP_BIN_NAME='TGLFEP_driver'
export TGLFEP_BIN_DIR=/global/homes/b/bassem/ips-fastran/examples/tglf-ep/
export TGLFEP_BIN_NAME=TGLFEP_driver

echo $TGLFEP_BIN_DIR
echo $TGLFEP_BIN_NAME

export SHOT_NUMBER=000001
export TIME_ID=00001

ips.py --simulation=fastran_scenario.config --platform=perlmutter_cpu_node.conf --log=ips.log 1> ips.out 2> ips.err &

conda deactivate
