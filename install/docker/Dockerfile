FROM base:latest

RUN apt-get install --yes  \
        libx11-dev           
RUN apt-get clean all

# ips framework
COPY IPS-framework /code/IPS-framework
RUN cd /code/IPS-framework && python setup.py install

COPY pyps /code/pyps
ENV PYTHONPATH=/code/pyps:$PYTHONPATH
RUN cd /code/pyps/test && python3 test1.py && python3 test2.py && python3 test3.py

# ips-fastran
COPY ips-fastran /code/ips-fastran
RUN cd /code/ips-fastran/src && python setup.py install
RUN cd /code/ips-fastran/lib && python setup.py install

# ips-eped
COPY ips-eped /code/ips-eped
RUN cd /code/ips-eped/src && python setup.py install

# copy physics codes
COPY binaries /code/binaries
COPY data /code/data

# environment variables

ENV PS_BACKEND=pyps

ENV DATA_DIR="/code/data"

ENV FASTRAN_BIN_DIR="/code/binaries"
ENV FASTRAN_BIN_NAME=xfastran_ver0.98
ENV FASTRAN_SERIAL_BIN_NAME=xfastran_ver0.98_ser

ENV EFIT_BIN_DIR="/code/binaries"
ENV EFIT_BIN_NAME=efitd90_static

ENV EFIT_TABLE_DIR="/code/data/greentable/"
ENV EFIT_INPUT_DIR="/code/data/greentable/"

ENV ESC_BIN_DIR="/code/binaries"
ENV ESC_BIN_NAME=xesc

ENV NUBEAM_BIN_DIR="/code/binaries"
ENV NUBEAM_BIN_NAME=nubeam_comp_exec

ENV NFREYA_BIN_DIR="/code/binaries"
ENV NFREYA_BIN_PATH="/code/binaries"
ENV NFREYA_BIN_NAME=onetwo_129_201
ENV NFREYA_DATA_ROOT="/code/data/onetwo"

ENV GENRAY_BIN_DIR="/code/binaries"
ENV GENRAY_BIN_NAME=xgenray.intel.edison

ENV TORAY_BIN_DIR="/code/binaries"
ENV TORAY_BIN_NAME=xtoray
