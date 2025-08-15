#!/bin/bash -l

#SBATCH --image=docker:registry.services.nersc.gov/rwp53/ips-massive-serial:dev
#SBATCH --volume="/global/cscratch1/sd/rwp/tmpfiles:/tmp:perNodeCache=size=20G"

#SBATCH -p debug
#SBATCH --nodes=2
#SBATCH -t 00:10:00
#SBATCH -C haswell

#SBATCH -J ips_fastran
#SBATCH -e ips.err
#SBATCH -o ips.out

export BIN_DIR=/global/common/software/atom/cori/binaries
module load python
conda activate massiveparallel

ips.py --config=ips_massive_serial_global_shifter.config --platform=cori_haswell.conf --log=ips.log
