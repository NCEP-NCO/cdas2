#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

f=r2_prepdate
target=cdas2_v1_prepdate

FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -assume byterecl

${target}:	$f.f
	${FC} $(FFLAGS) -o ${target} $f.f ${bufrlib}

clean:
	touch ${target} ; rm ${target}
