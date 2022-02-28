#LIBS=-L/nwprod/lib -lbacio_4
LIBS=${BACIO_LIB4}

FC=ifort
# FC=ftn
CC=icc
# CC=cc


#w3=getgbss.o rdgb.o baread.o baopen.o getgir.o ixgb.o skgb.o \
#w3fi01.o w3fi63.o w3fi83.o pdseup.o mova2i.o gbytes.o gbytes_char.o

w3=getgbss.o rdgb.o getgir.o ixgb.o skgb.o \
w3fi01.o w3fi63.o w3fi83.o pdseup.o mova2i.o gbytes.o gbytes_char.o

# w3=rdgb.o getgbss.o

# FFLAGS=-O2 -qrealsize=8 -qmaxmem=-1
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens

# w3flags=-O2 -qrealsize=4 -qmaxmem=-1
w3flags=-O2 -assume byterecl -assume protect_parens -convert big_endian

CFLAGS=-c -O2

target=cdas2_v1_sfc

f=r2_sfc

${target}:	$f.f $(w3)
	${FC} -o ${target} ${FFLAGS} $f.f ${w3} ${LIBS}

.f.o:
	${FC} ${w3flags} -c $<

.c.o:
	${CC} ${CFLAGS} $<

clean:
	touch ${target} junk.o ; rm ${target} *.o
