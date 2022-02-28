#w3lib=-L/nwprod/lib -lw3nco_d
w3lib=${W3NCO_LIBd}

#bufrlib=-L/nwprod/lib -lbufr_d_64
bufrlib=${BUFR_LIBd}

obj=chlml.o conmc.o dtast.o dvast.o fulldivt.o  g2s0.o gdcrdn.o gdcrdp.o \
	genqsat.o getbaln.o getlalo.o getpln.o  glbsoi2.o grad2s.o gtbdivt.o \
	gtbhalf.o hoper.o hopers.o htoper.o  inguess.o inguessv.o  initps2.o \
	initqpw2.o initsat6.o  initt3.o initw2.o intps2.o intqpw2.o  intrp2.o \
	intrp3.o intt3.o intw2.o  m1glat.o m1ipqr.o m1poly.o m1rcons2.o  pcgsoi5.o  \
	prepp.o preppw.o prepq.o preps.o prept.o prepw.o prnon85.o qoper.o \
	qtoper.o  rdfact.o rdgesc.o rdprep7.o rdtest9.o  residw.o respsf.o \
	respw.o resq.o ressat.o restmp.o  s2g0.o s2grad.o s2gvec.o s2mg2x.o \
	satc2.o satcov.o satop46.o  setuprhs5.o sprp2.o sprqpw2.o sprs6.o sprt2.o \
	spruv2.o  ssi2.o tg2s0.o tgrad2s.o ts2g0.o ts2grad.o ts2gvec.o wranlc.o  \
	fft99m.o sdot.o sgemv.o minv.o isrchfle.o isrchfge.o w3fs03.o w3fs11.o \
	w3fs21.o mova2i.o

FC=ifort
#FC=ftn
FFLAGS=-O2 -r8 -convert big_endian -assume byterecl -assume protect_parens

CC=icc
# CC=cc
CFLAGS=-O2


w3=iw3jdn.o

f=r2_ssi
target=cdas2_v1_ssi

${target}:	$(obj) $(w3)
	${FC} -o ${target} $(obj) $(bufrlib) $(w3)

.f.o:
	${FC} ${FFLAGS} -c $<

.c.o:
	${CC} ${CFLAGS} -c $<

clean:
	rm -f *.o ${target}
