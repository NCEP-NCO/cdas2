f=r2_sigdate
target=cdas2_v1_sigdate.gnu
FC=gfortran
FFLAGS=-O2 -fdefault-real-8 -fconvert=big-endian

${target}:	$f.f
	${FC} $(FFLAGS) -o ${target} $f.f
