SHELL=/bin/sh

#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

PROG=r2_fcst
TARGET = cdas2_v1_fcst
FFLAGS = -O2  -r8 -convert big_endian -assume byterecl -assume protect_parens -g -traceback

FC=ifort
# FC=ftn

${TARGET}:   ${PROG}.o 
	${FC} $(FFLAGS) -o ${TARGET} $(PROG).o $(w3lib)

${PROG}.o:   ${PROG}.f
	${FC} -c $(FFLAGS) ${PROG}.f

clean:
	touch ${TARGET} junk.o
	rm ${TARGET} *.o
