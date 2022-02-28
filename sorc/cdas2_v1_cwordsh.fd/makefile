# LIBS = -L/nwprod/lib -lw3nco_4 -lbufr_v10.2.0_4_64
LIBS=${BUFR_LIB4} ${W3NCO_LIB4}
FC=ifort

PROG=cdas2_v1_cwordsh

${PROG}:	cwordsh.f
	${FC} -O3 -assume noold_ldout_format -Xlinker -s -o ${PROG} cwordsh.f $(LIBS)

clean:
	rm ${PROG}
