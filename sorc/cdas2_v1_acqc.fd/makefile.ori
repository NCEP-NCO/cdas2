#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

FC=ifort
# FC=ftn
# FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens

f=r2_acqc
target=cdas2_v1_acqc

${target}:	$f.o
	${FC} $(FFLAGS) -o ${target} $f.o w3fi04.f ${w3lib} ${bufrlib}

.f.o:
	${FC} $(FFLAGS) -c $<
	
clean:
	touch junk.o ${target}; rm *.o ${target}
