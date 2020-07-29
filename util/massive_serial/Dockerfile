FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN \
    apt-get update        && \
    apt-get install --yes    \
        build-essential      \
        gfortran             \
        python3-dev          \
        wget                 \
        cmake                \
        doxygen              \
        git                  \
        python3-numpy        \
        python3-scipy        \
        python3-matplotlib   \
        python3-netcdf4      \
        libnetcdff-dev       \
        libopenblas-dev   && \
    apt-get clean all

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN mkdir /code

ARG GITHUB_USER
ARG GITHUB_TOKEN

RUN cd /code && \
    git clone --depth 1 https://github.com/HPC-SimTools/IPS-framework.git ips-framework              && \
    git clone --single-branch --branch massive-serial --depth 1 https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/ORNL-Fusion/ips-fastran.git

ENV IPS_DIR /code

# byte-compile everything we might need, exclude a few files that can't compile
RUN cd /code && py3compile ips-framework/framework/src           \
                           ips-framework/framework/utils/HTML.py \
                           ips-fastran/lib/python/Namelist.py    \
                           ips-fastran/src                       \
                           -X ".*/elite_io.py"                   \
                           -X ".*/modeleq_constraint_pscale.py"  \
                           -X ".*/pest3.py"

ENV IPS_PATH $IPS_DIR/ips-framework
ENV IPS_ROOT $IPS_PATH
ENV PATH $IPS_PATH/framework/src:/code/bin:$PATH
ENV PYTHONPATH $IPS_PATH/framework/src:$IPS_PATH/framework/utils:$PYTHONPATH

ENV FASTRAN_ROOT="${IPS_DIR}/ips-fastran"
ENV PYTHONPATH="$FASTRAN_ROOT/lib:$PYTHONPATH"
ENV PYTHONPATH="$FASTRAN_ROOT/lib/python:$PYTHONPATH"
ENV PYTHONPATH="$FASTRAN_ROOT/src:$PYTHONPATH"

RUN mkdir -p $FASTRAN_ROOT/lib/omfit/classes/ && touch $FASTRAN_ROOT/lib/omfit/classes/__init__.py
COPY fluxSurface.py $FASTRAN_ROOT/lib/omfit/classes/
RUN chmod 644 $FASTRAN_ROOT/lib/omfit/classes/fluxSurface.py
RUN python3 -c "from omfit.classes.fluxSurface import boundaryShape"

COPY pyps.tar.gz /code
RUN cd /code && tar xf pyps.tar.gz && chmod 755 pyps_new && chmod 644 pyps_new/*.py
ENV PYTHONPATH="${IPS_DIR}/pyps_new:$PYTHONPATH"
RUN cd /code/pyps_new/pyps_test3 && python3 test1.py && python3 test2.py && python3 test3.py

ENV PS_BACKEND=pyps

RUN mkdir /code/bin

COPY efitd90_static /code/bin
RUN chmod 755 /code/bin/efitd90_static
ENV EFIT_BIN_DIR="${IPS_DIR}/bin"
ENV EFIT_BIN_NAME=efitd90_static

ENV EFIT_TABLE_DIR="/global/common/software/atom/cori/binaries/efit/greentable/"
ENV EFIT_INPUT_DIR="/global/common/software/atom/cori/binaries/efit/greentable/"

COPY xfastran_ver0.93_ser /code/bin
RUN chmod 755 /code/bin/xfastran_ver0.93_ser
RUN cd /code/bin && ln -s xfastran_ver0.93_ser xfastran
ENV FASTRAN_BIN_DIR="${IPS_DIR}/bin"
ENV FASTRAN_BIN_NAME=xfastran
ENV FASTRAN_SERIAL_BIN_NAME=xfastran_ver0.93_ser

# https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
ENV OPENBLAS_NUM_THREADS=1

COPY ips_continuous_copy.sh /code/bin