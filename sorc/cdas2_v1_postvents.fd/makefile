#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}


f=r2_postvents
target=cdas2_v1_postvents

FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl

${target}:	$f.f
	${FC} -o ${target} ${FFLAGS} $f.f ${w3lib} ${bufrlib}

clean:
	touch ${target} ; rm ${target}
