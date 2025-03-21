#!/bin/bash -l

module load python
source activate /global/common/software/atom/perlmutter/cesol/conda/dev
export PYTHONPATH=$HOME/ips-fastran/src:$PYTHONPATH # <== set to your local path

export SHOT_NUMBER=000001
export TIME_ID=00001

ips.py --simulation=fastran_scenario.config --platform=$MACHINE_CONFIG_SERIAL --log=ips.log 1> ips.out 2> ips.err &

conda deactivate
