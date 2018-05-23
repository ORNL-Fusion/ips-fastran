LIBHOME = ${IPS_ROOT}/ntcc 

FC = gfortran
FFLAGS = -c -O
INCLUDE = -I. -I${LIBHOME}/mod #-I${NETCDFDIR}/include

NTCCLIB = -L${LIBHOME}/lib \
            -lplasma_state -lps_xplasma2 -lplasma_state_kernel \
            -lxplasma_debug -lxplasma2 -lgeqdsk_mds -lmdstransp \
            -lnscrunch -lfluxav -lr8bloat -ltrgraf -lureadsub -lvaxonly \
            -lpspline -lezcdf -llsode -llsode_linpack -lelvislib \
            -lsg -ljc -lsmlib -lcomput -lportlib  -L/opt/local/lib -lnetcdf -llapack -lblas


NTCC = ${NTCCLIB}

CFLAGS_ALL = ${CFLAGS}
CDEFS_ALL  = ${CDEFS}
FFLAGS_ALL = ${FFLAGS} ${INCLUDE}
FDEFS_ALL  = ${FDEFS}

SRC_C   = 
SRC_F   = 
SRC_F90 = pstool.f90

OBJS = ${SRC_C:.c=.o} ${SRC_F90:.f90=.o} ${SRC_F:.f=.o}

%.o : %.c 
	${CC} ${CFLAGS_ALL} ${CDEFS_ALL} $*.c

%.o : %.f
	${FC} ${FFLAGS_ALL} ${FDEFS_ALL} $*.f

%.o : %.f90
	${FC} ${FFLAGS_ALL} ${FDEFS_ALL} $*.f90

pstool : ${OBJS}
	${FC} -o pstool ${OBJS} ${NTCC}

clean:
	rm *.o
	rm *.a
	rm *.mod

