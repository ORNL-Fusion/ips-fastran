SWIG = swig
SWIG_FLAGS = -python -c++
CXX_GCC = g++
CXX = icc
CFLAGS = -g -c -fPIC 

ATOM = /project/projectdirs/atom/atom-install-edison
LOCAL = ${ATOM}/cesol

NTCCDIR =  ${LOCAL}/local/ntcc_intel
NETCDFDIR = ${NETCDF_DIR}/lib
LDFLAGS = -shared -fPIC
LIBS = -L. \
    -L${NTCCDIR}/lib \
    -L${NETCDFDIR} \
    -lccps_lib -lplasma_state -lps_xplasma2 -lplasma_state_kernel -lxplasma2 \
    -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch -lfluxav -lr8bloat -lpspline -lezcdf \
    -llsode -llsode_linpack -lsmlib -lcomput -lportlib \
    -lnetcdff -lnetcdf \
    -mkl -lifcore -lc -L/usr/common/software/python/3.7-anaconda-2019.10/lib -lpython3.7m

INCDIR = -I${NTCCDIR}/include_ntcc -I/usr/common/software/python/3.7-anaconda-2019.10/include/python3.7m


all: _pyplasmastate.so

pyplasmastate_wrap.cxx: pyplasmastate.i 
	${SWIG} ${SWIG_FLAGS} $? 

pyplasmastate_wrap.o: pyplasmastate_wrap.cxx 
	${CXX} ${CFLAGS} ${INCDIR} $?

_pyplasmastate.so: pyplasmastate_wrap.o
	${CXX} ${LDFLAGS} $? ${LIBS} -o $@

clean:
	rm -f *.o
	rm -f *.pyc
	rm -f _pyplasmastate.so
	rm -f pyplasmastate_wrap.cxx
	rm -f pyplasmastate.py
