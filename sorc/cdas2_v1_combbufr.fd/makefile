#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl

f=r2_combbufr
target=cdas2_v1_combbufr

${target}:	$f.f
	${FC} $(FFLAGS) -o ${target} $f.f ${bufrlib}

clean:
	touch ${target} ; rm ${target}
