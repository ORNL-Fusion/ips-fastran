FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN \
    apt-get update        && \
    apt-get install --yes    \
        build-essential      \
        gfortran             \
        wget                 \
        git                  \
        libnetcdff-dev       \
        libopenblas-dev      \
        &&                   \
    apt-get clean all

ARG mpich=3.3
ARG mpich_prefix=mpich-$mpich

RUN \
    wget https://www.mpich.org/static/downloads/$mpich/$mpich_prefix.tar.gz && \
    tar xvzf $mpich_prefix.tar.gz                                           && \
    cd $mpich_prefix                                                        && \
    ./configure                                                             && \
    make -j 4                                                               && \
    make install                                                            && \
    make clean                                                              && \
    cd ..                                                                   && \
    rm -rf $mpich_prefix

RUN /sbin/ldconfig

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

ENV PATH="/opt/conda/bin:${PATH}"
ARG PATH="/opt/conda/bin:${PATH}"

RUN conda install python=3.8 numpy scipy matplotlib netcdf4 dask

RUN python3 -m pip install mpi4py

RUN mkdir /code

RUN python -m pip install ipsframework

COPY greentable.tar.gz /code
RUN cd /code && tar xf greentable.tar.gz && rm greentable.tar.gz && chmod 755 greentable && chmod 644 greentable/*

ARG GITHUB_USER
ARG GITHUB_TOKEN

RUN cd /code && git clone --single-branch --branch update-massive-serial --depth 1 https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/ORNL-Fusion/ips-fastran.git
RUN cd /code/ips-fastran/src && python -m pip install .
RUN cd /code/ips-fastran/lib/namelist && python -m pip install .
RUN rm -r /code/ips-fastran

RUN mkdir -p /opt/conda/lib/python3.8/site-packages/omfit_classes/ && touch /opt/conda/lib/python3.8/site-packages/omfit_classes/__init__.py
COPY fluxSurface.py /opt/conda/lib/python3.8/site-packages/omfit_classes/
RUN chmod 644 /opt/conda/lib/python3.8/site-packages/omfit_classes/fluxSurface.py
RUN python -c "from omfit_classes.fluxSurface import boundaryShape"

COPY pyps_38.tar.gz /code
RUN cd /code && tar xf pyps_38.tar.gz && chmod 755 pyps && chmod 644 pyps/*.py && rm pyps_38.tar.gz
ENV PYTHONPATH="/code/pyps:$PYTHONPATH"
RUN cd /code/pyps/test && python3 test1.py && python3 test2.py && python3 test3.py

ENV PS_BACKEND=pyps

RUN mkdir /code/bin

ENV BIN_DIR=/code/bin
ENV PATH /code/bin:$PATH

COPY efitd90_static /code/bin
RUN chmod 755 /code/bin/efitd90_static
ENV EFIT_BIN_DIR=/code/bin
ENV EFIT_BIN_NAME=efitd90_static

ENV EFIT_TABLE_DIR=/code/greentable/
ENV EFIT_INPUT_DIR=/code/greentable/

COPY xfastran_ver0.93_ser /code/bin
RUN chmod 755 /code/bin/xfastran_ver0.93_ser
RUN cd /code/bin && ln -s xfastran_ver0.93_ser xfastran
ENV FASTRAN_BIN_DIR=/code/bin
ENV FASTRAN_BIN_NAME=xfastran
ENV FASTRAN_SERIAL_BIN_NAME=xfastran_ver0.93_ser

# https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
ENV OPENBLAS_NUM_THREADS=1

COPY ips_continuous_copy.sh /code/bin
