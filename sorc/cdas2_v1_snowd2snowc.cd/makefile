#LIB=-L/nwprod/lib/ -lgribw
#INCLUDE=-I/nwprod/lib/sorc/gribw/include
#LIB=${GRIBW_LIB} -I

GRIBW_SRC=../gribw_lib.cd
LIB=-L${GRIBW_SRC} -lgribw
INCLUDE=-I${GRIBW_SRC}

CC=icc

CMD=cdas2_v1_snowd2snowc

${CMD}:	snowd2snowc.c
	${CC} -o ${CMD} -O2 ${INCLUDE} snowd2snowc.c ${LIB} -lm

clean:
	rm -f *.o ${CMD}
