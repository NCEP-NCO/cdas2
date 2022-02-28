FC=ifort
# FC=ftn

f=r2_sigdate
target=cdas2_v1_sigdate
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl

${target}:	$f.f
	${FC} $(FFLAGS) -o ${target} $f.f

clean:
	touch ${target}
	rm ${target}
