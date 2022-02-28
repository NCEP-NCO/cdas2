SRCS= gribify_rt_p_pingrain.c
CMD= cdas2_v1_gribpentad
INCLUDE=gribwlib.h pdstool.h gdstool.h gribw_time.h 
#LIBS= -L/nwprod/lib -lgribw
# LIBS=${GRIBW_LIB}
LIBS= -L../gribw_lib.cd -lgribw

CC=icc
# CC=cc
CFLAGS=-O2

$(CMD):	$(SRCS) $(INCLUDE)
	$(CC) -o $(@) ${CFLAGS}  $(SRCS) $(LIBS) -lm

clean:
	touch ${CMD} ; rm ${CMD}
