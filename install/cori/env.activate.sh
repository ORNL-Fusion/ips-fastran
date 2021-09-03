#!/bin/sh

export BIN_DIR=/global/common/software/atom/cori/binaries
export DATA_DIR=/global/common/software/atom/cori/data
export EFIT_TABLE_DIR=$DATA_DIR/greentable/
export EFIT_INPUT_DIR=$DATA_DIR/greentable/
export WGEQDSK_BIN_DIR=$BIN_DIR/wgeqdsk/default
export WGEQDSK_BIN_NAME=wgeqdsk
export FASTRAN_BIN_DIR=$BIN_DIR/fastran/default
export FASTRAN_BIN_NAME=xfastran_ver0.97
export FASTRAN_SERIAL_BIN_NAME=xfastran_ver0.93_ser
export EFIT_BIN_DIR=$BIN_DIR/efit/static
export EFIT_BIN_NAME=efitd90_static
export ESC_BIN_DIR=$BIN_DIR/esc/default
export ESC_BIN_NAME=xesc
export NUBEAM_BIN_DIR=$BIN_DIR/nubeam/default
export NUBEAM_BIN_NAME=mpi_nubeam_comp_exec
export TORAY_BIN_DIR=$BIN_DIR/toray/default
export TORAY_BIN_NAME=xtoray
export CURRAY_BIN_DIR=$BIN_DIR/curray/default
export CURRAY_BIN_NAME=xcurray
export GENRAY_BIN_DIR=$BIN_DIR/genray/default
export GENRAY_BIN_NAME=xgenray.intel.edison
export NFREYA_BIN_PATH=$BIN_DIR/onetwo/default
export NFREYA_BIN_NAME=onetwo_129_201
export NFREYA_DATA_ROOT=$DATA_DIR/onetwo
export PSTOOL=$BIN_DIR/$PSTOOL_BIN_NAME
export PS_BACKEND=pyps

export DAKOTA_ROOT=/global/common/software/atom/cori/dakota
export PATH=$DAKOTA_ROOT/bin:$PATH
