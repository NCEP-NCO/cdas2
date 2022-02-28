#LIBS=-L/nwprod/lib -lw3nco_d -lw3emc_d -lbufr_d_64
LIBS=${W3NCO_LIBd} ${W3EMC_LIBd} ${BUFR_LIBd}

f=r2_prevents
target=cdas2_v1_prevents

FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl

${target}:	$f.f
	${FC} -o ${target} ${FFLAGS} $f.f ${LIBS} 

clean:
	touch ${target} ; rm ${target}
