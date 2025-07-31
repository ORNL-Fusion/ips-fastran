echo $CONDA_PREFIX
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
touch $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
touch $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
cp $MACHINE/env.activate.sh $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
cp $MACHINE/env.deactivate.sh $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

