#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

f=r2_precipadj
target=cdas2_v1_precipadj
FC=ifort
# FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl

${target}:	$f.f
	${FC} $(FFLAGS) -o ${target} $f.f ${w3lib}

clean:
	touch ${target} ; rm ${target}

