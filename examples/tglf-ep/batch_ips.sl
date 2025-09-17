#!/bin/bash

#SBATCH --qos=regular
#SBATCH -A m808
#SBATCH --constraint=cpu
#SBATCH -o ./run.out.ips
#SBATCH -e ./run.err.ips
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH -J ips_tglfep

#!/bin/bash -l
source /usr/share/lmod/lmod/init/bash
           
module load python
source activate /global/common/software/atom/perlmutter/cesol/conda/dev
export PYTHONPATH=/global/homes/b/bassem/ips-fastran/src:$PYTHONPATH
echo $PYTHONPATH

export EP_BIN_DIR='/global/common/software/atom/gacode_add_TGLFEP'
export TGLFEP_BIN_NAME='TGLF-EP/TGLFEP_driver'
export ALPHA_BIN_NAME='Alpha/Alpha_driver'
#export TGLFEP_BIN_DIR=$SCRATCH/tglf-ep/
#export TGLFEP_BIN_NAME=TGLFEP_driver

echo $TGLFEP_BIN_DIR
echo $TGLFEP_BIN_NAME

export SHOT_NUMBER=000001
export TIME_ID=00001

ips.py --simulation=fastran_scenario.config --platform=$MACHINE_CONFIG --log=ips.log

conda deactivate


