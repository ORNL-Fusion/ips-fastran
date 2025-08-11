#!/bin/bash

#SBATCH --qos=debug
#SBATCH -A m808
#SBATCH --constraint=cpu
#SBATCH -o ./run.out.ips
#SBATCH -e ./run.err.ips
#SBATCH --ntasks-per-node=125
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00
#SBATCH -J ips_tglfep

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

ips.py --simulation=fastran_scenario.config --platform=$MACHINE_CONFIG --log=ips.log 1> ips.out 2> ips.err &

conda deactivate


