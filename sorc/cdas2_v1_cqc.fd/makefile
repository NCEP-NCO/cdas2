#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl  -traceback -save -i8

f=r2_cqc
target=cdas2_v1_cqc

${target}:	$f.F
	${FC} $(FFLAGS) -o ${target} $f.F ${w3lib} ${bufrlib}

clean:
	touch ${target} ; rm ${target}
