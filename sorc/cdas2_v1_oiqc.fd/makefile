SHELL=/bin/sh

#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens -traceback -qopenmp
# debug FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens -traceback -check all
FC=ifort
# FC=ftn

PROG=r2_oiqc
TARGET = cdas2_v1_oiqc
OBJS = r2_oiqc.o ufbrew.o

${TARGET}:   ${OBJS} 
	${FC} $(FFLAGS) -o ${TARGET} $(OBJS) $(w3lib) $(bufrlib)

${PROG}.o: ${PROG}.F
	${FC} $(FFLAGS) -c  ${PROG}.F
	
ufbrew.o: ufbrew.f
	${FC} $(FFLAGS2) -c  ufbrew.f

clean:
	touch junk.o ${TARGET} ; rm *.o ${TARGET}
