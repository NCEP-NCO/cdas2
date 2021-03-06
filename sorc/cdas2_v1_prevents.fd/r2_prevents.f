      SUBROUTINE FAX(IFAX,N,MODE)
       SAVE
      DIMENSION IFAX(10)
C
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
C     IFAX(1) CONTAINS NUMBER OF FACTORS
C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
      SUBROUTINE FFTRIG(TRIGS,N,MODE)
       SAVE
      DIMENSION TRIGS(1)
C
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5   E   0*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5  E  0*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5  E  0*FLOAT(I-1)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5  E  0*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=2.0  E  0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5  E  0*DEL
      DO 50 I=2,N
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
       SAVE
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,
     *     SIN60/0.866025403784437/
C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
C     CODING FOR FACTOR 2
C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 3
C
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                       -(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                       +(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                       +(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                       -(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                      -(SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                      +(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     *    S1*((A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                      -(SIN60*(B(IB+I)-B(IC+I))))
     *   +C1*((B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                      +(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                      +(SIN60*(B(IB+I)-B(IC+I))))
     *   -S2*((B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                      -(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-0.5  E  0*(A(IB+I)+A(IC+I)))
     X                                      +(SIN60*(B(IB+I)-B(IC+I))))
     *   +C2*((B(IA+I)-0.5  E  0*(B(IB+I)+B(IC+I)))
     X                                      -(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 4
C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 5
C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99M(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
       SAVE
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
      NFAX=IFAX(1)
      NX=N
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN ZEROS AT END
  110 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
       SAVE
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C     (SPECTRAL TO GRIDPOINT TRANSFORM)
C
      DIMENSION A(N),WORK(N),TRIGS(N)
C
      NH=N/2
      NX=N
      INK=INC+INC
C
C     A(0)   A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
      DO 10 L=1,LOT
      WORK(JA)=A(IA)
      WORK(JB)=A(IA)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *    (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *    (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      DO 40 L=1,LOT
      WORK(JA)=2.0  E  0*A(IA)
      WORK(JA+1)=-2.0  E  0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
       SAVE
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C     (GRIDPOINT TO SPECTRAL TRANSFORM)
C
      DIMENSION WORK(N),A(N),TRIGS(N)
C
      NH=N/2
      NX=N
      INK=INC+INC
C
C     A(0)   A(N/2)
      SCALE=1.0  E  0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JA+INC)=0.0  E  0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      SCALE=0.5  E  0*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0  E  0*SCALE
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   .===================================================================
C   | DESCRIPTION | NAME   | BITS | MAX VALUE  | DIMENSIONS
C   |=============*========*======*============*========================
C   | ZS          | ZS     | WORD | REAL       |  256 , 129
C   |-------------+--------+------+------------+------------------------
C   | PS          | PS     | WORD | REAL       |  256 , 129
C   |-------------+--------+------+------------+------------------------
C   | T           | T      | WORD | REAL       |  256 , 129 , 28
C   |-------------+--------+------+------------+------------------------
C   | U           | U      | WORD | REAL       |  256 , 129 , 28
C   |-------------+--------+------+------------+------------------------
C   | V           | V      | WORD | REAL       |  256 , 129 , 28
C   |-------------+--------+------+------------+------------------------
C   | Q           | Q      | WORD | REAL       |  256 , 129 , 28
C   `=============^========^======^============^========================
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION ZS(I,J)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      COMMON /STOR1/IAR(IM,JM)
      REAL IAR,ZS
      ZS = IAR(I,J)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE ZSP(I,J,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      COMMON /STOR1/IAR(IM,JM)
      REAL IAR,V
      IAR(I,J) = V
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION PS(I,J)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      COMMON /STOR2/IAR(IM,JM)
      REAL IAR,PS
      PS = IAR(I,J)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PSP(I,J,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      COMMON /STOR2/IAR(IM,JM)
      REAL IAR,V
      IAR(I,J) = V
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION T(I,J,K)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR3/IAR(IM,JM,KM)
      REAL IAR,T
      T = IAR(I,J,K)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TP(I,J,K,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR3/IAR(IM,JM,KM)
      REAL IAR,V
      IAR(I,J,K) = V
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION U(I,J,K)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR4/IAR(IM,JM,KM)
      REAL IAR,U
      U = IAR(I,J,K)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UP(I,J,K,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR4/IAR(IM,JM,KM)
      REAL IAR,V
      IAR(I,J,K) = V
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION V(I,J,K)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR5/IAR(IM,JM,KM)
      REAL IAR,V
      V = IAR(I,J,K)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE VP(I,J,K,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR5/IAR(IM,JM,KM)
      REAL IAR,V
      IAR(I,J,K) = V
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION Q(I,J,K)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR6/IAR(IM,JM,KM)
      REAL IAR,Q
      Q = IAR(I,J,K)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE QP(I,J,K,V)
      PARAMETER (IM =  256 )
      PARAMETER (JM =  129 )
      PARAMETER (KM =  28 )
      COMMON /STOR6/IAR(IM,JM,KM)
      REAL IAR,V
      IAR(I,J,K) = V
      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE ALLFAC(PLN,DEP,INX,FAC,SNL,CSL,RSC)
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
      COMMON /GESPRM/ JCAP,JCAP1,JCAP2,JCAP1X2,MDIMA,MDIMB,MDIMC
 
C-CRA DIMENSION PLN(JCAP1,2,JMAX),DEP(MDIMC,2),INX(MDIMC,2)
C-CRA DIMENSION FAC(0:JCAP1,0:JCAP,3),SNL(IMAX),CSL(IMAX),RSC(JMAX)
 
      DIMENSION PLN( 63 ,2, 129 )
      DIMENSION DEP( 4158 ,2),INX( 4158 ,2)
      DIMENSION FAC(0: 63 ,0: 62 ,3)
      DIMENSION SNL( 256 ),CSL( 256 ),RSC( 129 )
 
      DATA PI180/.0174532/
 
C----------------------------------------------------------------------
      FA(X,Y)   = 6.3712E6 * SQRT((X**2-Y**2)/(4*X**2-1))/X
      FB(X,Y)   = 6.3712E6 * Y/(X**2+X)
      XLAT(J)   = (J-1)*DLAT-90.
      COLR(J)   = (90.0-ABS(XLAT(J)))*.0174532
      DEPS(X,Y) = SQRT((X**2-Y**2)/(4.0*X**2-1.0))
C-----------------------------------------------------------------------
 
C  CLEAR THE ARRAYS
C  ----------------
 
      DO I=1,JCAP1*2*JMAX
        PLN(I,1,1) = 0
      ENDDO
      DO I=1,MDIMC*2
        DEP(I,1) = 0
        INX(I,1) = 0
      ENDDO
      DO K=1,3
        DO J=0,JCAP
          DO I=0,JCAP1
            FAC(I,J,K) = 0
          ENDDO
        ENDDO
      ENDDO
      DO I=1,IMAX
        SNL(I) = 0
        CSL(I) = 0
      ENDDO
      DO I=1,JMAX
        RSC(I) = 0
      ENDDO
 
C  COMPUTE THE TRANSPOSE INDEXES FOR MDIMA AND MDIMC
C  -------------------------------------------------
 
      L = 1
      DO M=1,JCAP1
      DO N=0,JCAP1-M
      IND = N*(JCAP1X2-N+1) + 2*M - 1
      INX(L  ,1) = IND
      INX(L+1,1) = IND+1
      L=L+2
      ENDDO
      ENDDO
 
      L = 1
      DO M=1,JCAP1
      DO N=0,JCAP1-M+1
      IF(N.EQ.0) IND = 2*M-1
      IF(N.EQ.1) IND = 2*(JCAP1+M)-1
      IF(N.GT.1) IND = N*(JCAP1X2-N+3)+2*M-3
      INX(L  ,2) = IND
      INX(L+1,2) = IND+1
      L=L+2
      ENDDO
      ENDDO
 
C  COMPUTE THE PLN FACTORS FOR EACH LATITUDE EXCEPT THE POLES
C  ----------------------------------------------------------
 
      DO J=2,JMAX-1
      SINLAT = COS(COLR(J))
      COS2   = 1.0-SINLAT**2
      PROD   = 1.0
      DO N=1,JCAP1
      X = 2*N+1
      SRHP = SQRT(PROD*.5)
      PLN(N,1,J) = SRHP
      PLN(N,2,J) = SRHP*SINLAT*SQRT(X)
      PROD  = PROD*COS2*(X/(X-1.))
      ENDDO
      ENDDO
 
C  COMPUTE DEPS FOR MDIMA AND MDIMC
C  --------------------------------
 
      IATA = 1
      IATC = 1
      LEN  = JCAP
      DO N=0,JCAP1
      DO M=0,LEN
      DEPX = DEPS(FLOAT(N+M),FLOAT(M))
      DEP(IATA+2*M  ,1) = DEPX
      DEP(IATC+2*M  ,2) = DEPX
      DEP(IATA+2*M+1,1) = DEPX
      DEP(IATC+2*M+1,2) = DEPX
      ENDDO
      IATA = IATA+2*(JCAP1-N)
      IATC = IATC+2*(LEN+1)
      LEN  = LEN-MIN(N,1)
      ENDDO
 
C  THE DZTOUV FACTORS
C  ------------------
 
      DO M=0,JCAP
      DO N=M,JCAP1
      X = N
      Y = M
      IF(N.LT.JCAP              ) FAC(N,M,1) = FA(X+1,Y)
      IF(N.GT.0                 ) FAC(N,M,2) = FA(X,Y)
      IF(N.GT.0 .AND. N.LT.JCAP1) FAC(N,M,3) = FB(X,Y)
      ENDDO
      ENDDO
 
C  POLAR WIND FACTORS
C  ------------------
 
      DO I=1,IMAX
      ANG = (I-1)*DLON*PI180
      SNL(I) = SIN(ANG)
      CSL(I) = COS(ANG)
      ENDDO
 
C  RECIPROCALS OF SINES OF COLATITUDES
C  -----------------------------------
 
      DO J=2,JMAX-1
      RSC(J) = 1./SIN(COLR(J))
      ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE COF2GRD(LUN,IRET)
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/  DLAT,DLON,SL(100),SI(101)
      COMMON /GESPRM/ JCAP,JCAP1,JCAP2,JCAP1X2,MDIMA,MDIMB,MDIMC
 
C-CRA DIMENSION COF(MDIMC,2),GRI(MDIMC,2),GRD(IMAX,JMAX,2)
C-CRA DIMENSION FAC(0:JCAP1,0:JCAP,3),SNL(IMAX),CSL(IMAX),RSC(JMAX)
C-CRA DIMENSION INX(MDIMC,2),PLN(JCAP1,2,JMAX),DEP(MDIMC,2)
C-CRA DIMENSION IFAX(20),TRIGS(IMAX,2),WORK(IMAX,4)
 
      DIMENSION COF( 4158 ,2),GRI( 4158 ,2),GRD( 256 , 129 ,2)
      DIMENSION FAC(0: 63 ,0: 62 ,3),SNL( 256 )
      DIMENSION CSL( 256 ),RSC( 129 )
      DIMENSION INX( 4158 ,2),PLN( 63 ,2, 129 ),DEP( 4158 ,2)
      DIMENSION IFAX(20),TRIGS( 256 ,2),WORK( 256 ,4)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  INITIALIZE SOME STUFF
C  ---------------------
 
      CALL ALLFAC(PLN,DEP,INX,FAC,SNL,CSL,RSC)
C-CRA CALL FFTFAX(IMAX,IFAX,TRIGS)
      CALL    FAX(IFAX, IMAX,3)
      CALL FFTRIG(TRIGS,IMAX,3)
      IF(IFAX(1).EQ.0 .OR. IFAX(1).EQ.-99) GOTO 902
 
C  READ AND TRANSFORM FIELDS IN THE SPECTRAL FILE
C  ----------------------------------------------
 
      DO I=1,5
      NRD = 1
      LEV = KMAX
      NCF = MDIMA
      IF(I.LE.2) LEV = 1
      IF(I.EQ.4) NRD = 2
      IF(I.EQ.4) NCF = MDIMC
      DO L=1,LEV
      DO K=1,NRD
      READ(LUN,END=900,ERR=901) (COF(II,K),II=1,MDIMA)
      ENDDO
 
      IF(I.EQ.4) CALL DZTOUV(COF,FAC)
 
      DO J=2,JMAX-1
      CALL PLNSUM(COF,NCF,J,INX(1,NRD),PLN(1,1,J),DEP(1,NRD),GRI)
      DO K=1,NRD
C-CRA CALL RFFTMLT(GRI(1,K),WORK,TRIGS,IFAX,1,IMAX,IMAX,1,1)
      CALL FFT99M (GRI(1,K),WORK,TRIGS,IFAX,1,IMAX,IMAX,1,1)
      DO IM=1,IMAX
      GRD(IM,J,K) = GRI(IM,K)
      ENDDO
      ENDDO
      ENDDO
      CALL POLES(GRD,I,SNL,CSL,RSC)
      CALL GUSER(GRD,I,L)
      PRINT *,'I=',I,' L=',L
      CALL MAXMIN(GRD,IMAX,JMAX,IMAX,JMAX,NRD)
      ENDDO
      ENDDO
 
      RETURN
900   PRINT*,'COF2GRD - EOF READING GUESS  '
      IRET = -1
      RETURN
901   PRINT*,'COF2GRD - ERROR READING GUESS'
      IRET = -1
      RETURN
902   PRINT*,'COF2GRD - IDIM NOT FACTORABLE'
      IRET = -1
      RETURN
      END
C----------------------------------------------------------------------
C  THIS ROUTINE PERFORMS CONVERSION OF TRUE SCALAR RELATIVE
C  VORTICITY (Z), DIVERGENCE (D) TO PSEUDO SCALAR U AND V
C  IN SPECTRAL SPACE.
C
C  CALL DZTOUV(NLV,U,V)
C
C  NLV ......... NUMBER OF LEVELS TO CONVERT
C  U (INPUT) ... DIVERGENCE
C  V (INPUT) ... RELATIVE VORTICITY
C  U (OUTPUT) .. OUTPUT PSEUDO SCALAR ZONAL WIND (=UCOS(PHI))
C  V (OUTPUT) .. OUTPUT PSEUDO SCALAR MERID WIND (=VCOS(PHI))
C
C  HUA-LU PAN   27 FEBRUARY 1989
C  J WOOLLEN    14 MARCH    1993
C
C   MODIFIED TO RUN ON THE CRAY
C----------------------------------------------------------------------
      SUBROUTINE DZTOUV(COF,F)
 
      COMMON /GESPRM/ JCAP,JCAP1,JCAP2,JCAP1X2,MDIMA,MDIMB,MDIMC
 
C-CRA DIMENSION COF(MDIMC,2),F(0:JCAP1,0:JCAP,3)
C-CRA DIMENSION D(-1:MDIMA+4),Z(-1:MDIMA+4)
C-CRA DIMENSION U(MDIMC),V(MDIMC)
 
      DIMENSION COF( 4158 ,2),F(0: 63 ,0: 62 ,3)
      DIMENSION D(-1: 4032 +4),Z(-1: 4032 +4)
      DIMENSION U( 4158 ),V( 4158 )
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C  CONVERT D AND Z TO SCALED U AND V
C  ---------------------------------
 
      DO I=-1,MDIMA+4
        D(I) = 0
        Z(I) = 0
      ENDDO
 
      DO I=1,MDIMA
      D(I) = COF(I,1)
      Z(I) = COF(I,2)
      ENDDO
 
      MDZ = 0
      MUV = 0
 
      DO M=0,JCAP
      DO N=M,JCAP1
      MUV = MUV+2
      MDZ = MDZ+2
      U(MUV-1) =  Z(MDZ+1)*F(N,M,1)-Z(MDZ-3)*F(N,M,2)+D(MDZ  )*F(N,M,3)
      U(MUV  ) =  Z(MDZ+2)*F(N,M,1)-Z(MDZ-2)*F(N,M,2)-D(MDZ-1)*F(N,M,3)
      V(MUV-1) = -D(MDZ+1)*F(N,M,1)+D(MDZ-3)*F(N,M,2)+Z(MDZ  )*F(N,M,3)
      V(MUV  ) = -D(MDZ+2)*F(N,M,1)+D(MDZ-2)*F(N,M,2)-Z(MDZ-1)*F(N,M,3)
      ENDDO
      MDZ = MDZ-2
      ENDDO
 
      DO I=1,MDIMC
      COF(I,1) = U(I)
      COF(I,2) = V(I)
      ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE GESRES(LUN,IRET)
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
      COMMON /GESPRM/ JCAP,JCAP1,JCAP2,JCAP1X2,MDIMA,MDIMB,MDIMC
      COMMON /GESDAT/ FCLABL(4),IDATE,VDATE
      CHARACTER*8 FCLABL,IDATE,VDATE
 
      DIMENSION HEADR2(207)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      IRET = 0
 
C  READ SIGGES HEADERS
C  -------------------
 
C-MK  REWIND LUN
C-MK  READ(LUN,END=900,ERR=901) FCLABL
C-MK  READ(LUN,END=900,ERR=5  ) FHR,IH,IM,ID,IY,HEADR2
C-MK  goto 6
C-MK
C-MK5     continue
C-MK  print*,' attempting to read ges resolution from unit 5'
C-MK  read(5,*,end=901,err=901) jcap,kmax
C-MK
C-MK
 
      KMAX= 28
 
      REWIND LUN
      READ(LUN,END=900,ERR=901) FCLABL
      READ(LUN,END=900,ERR=901) FHR,IH,IM,ID,IY,(headr2(i),i=1,2*kmax+1)
C-MK
C-MK  headr2(202) = jcap
C-MK  headr2(203) = kmax
C-MK6     continue
 
C  EXTRACT HEADER INFO
C  -------------------
 
C-MK  JCAP  = HEADR2(202)
C-MK  KMAX  = HEADR2(203)
 
      JCAP= 62
 
      IF(KMAX.GT.100) GOTO 902
 
      DO L=1,KMAX
      SI(L) = HEADR2(L)
      SL(L) = HEADR2(KMAX+1+L)
      ENDDO
 
      SI(KMAX+1) = HEADR2(KMAX+1)
 
      PRINT *,'SI=',(SI(L),L=1,KMAX+1)
      PRINT *,'SL=',(SL(L),L=1,KMAX)
 
      JHR = FHR
      CALL W3FS03(IDI,IH,IY,IM,ID,0)
      CALL W3FS15(IDI,JHR,IDV)
      CALL W3FS03(IDV,JH,JY,JM,JD,1)
      WRITE(IDATE,'(4I2)') IY,IM,ID,IH
      WRITE(VDATE,'(4I2)') JY,JM,JD,JH
      DO I=1,8
      IF(IDATE(I:I).EQ.' ') IDATE(I:I) = '0'
      IF(VDATE(I:I).EQ.' ') VDATE(I:I) = '0'
      ENDDO
      WRITE(6,1) JHR,IDATE,VDATE
1     FORMAT(1X,'A ',I3,' HOUR FORECAST FROM ',A8,' VALID AT ',A8)
 
C  DEFINE THE OTHER RESOLUTION PARAMETERS
C  --------------------------------------
 
      JCAP1   = JCAP+1
      JCAP2   = JCAP+2
      JCAP1X2 = JCAP1*2
      MDIMA   = JCAP1*JCAP2
      MDIMB   = MDIMA/2+JCAP1
      MDIMC   = MDIMB*2
C
      IMAX    =  256
      JMAX    = IMAX/2+1
C
C EQUIVALENT TO STANDARD PLI PARAMETERS
C
C KMAX=LEVS
C IMAX=LONSSI
C JMAX=LATSSI
C JCAP1=JCAP1
C JCAP2=JCAP2
C JCAP1X2=TWOJ1
C MDIMA=LNUV
C MDIMB=LNUT
C MDIMC=LNUT2
C
      IF(IMAX.LT.JCAP1X2) GOTO 903
 
      DLAT  = 180./(JMAX-1)
      DLON  = 360./IMAX
 
      WRITE(6,2) JCAP,KMAX,DLAT,DLON
2     FORMAT(1X,'T',I3,' ',I2,' LEVELS -------> ',F3.1,' X ',F3.1)
 
      CALL COF2GRD(LUN,IRET)
 
      RETURN
900   PRINT*,'GESRES - EOF   READING GUESS ON LUN ',LUN
      IRET = -1
      RETURN
901   PRINT*,'GESRES - ERROR READING GUESS ON LUN ',LUN
      IRET = -1
      RETURN
902   PRINT*,'GESRES - KMAX TOO BIG = ',KMAX
      IRET = -1
      RETURN
903   PRINT*,'GESRES - IMAX TOO SMALL = ',IMAX
      IRET = -1
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLNSUM(COF,NCF,J,INX,PLN,DEP,GRD)
 
      COMMON /GESPRM/ JCAP,JCAP1,JCAP2,JCAP1X2,MDIMA,MDIMB,MDIMC
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
 
C-CRA DIMENSION COF(MDIMC,2),GRD(MDIMC,2)
C-CRA DIMENSION INX(MDIMC),PLN(JCAP1,2),DEP(MDIMC)
C-CRA DIMENSION QLN(MDIMC)
 
      DIMENSION COF( 4158 ,2),GRD( 4158 ,2)
      DIMENSION INX( 4158 ),PLN( 63 ,2),DEP( 4158 )
      DIMENSION QLN( 4158 )
 
C-----------------------------------------------------------------------
      XLAT(J)   = (J-1)*DLAT-90.
      COLR(J)   = (90.0-ABS(XLAT(J)))*.0174532
C-----------------------------------------------------------------------
 
      DO I=1,MDIMC*2
        GRD(I,1) = 0
      ENDDO
 
C  SET THE FIELD RELATED PARAMETERS
C  --------------------------------
 
      IF(NCF.EQ.MDIMA) THEN
         NCR = 1
         N0  = 1
      ELSE IF(NCF.EQ.MDIMC) THEN
         NCR = 2
         N0  = 0
      ELSE
         CALL SABORT('PLNSUM - UNKNOWN VALUE OF NCF')
      ENDIF
 
C  TRANSPOSE THE INPUT COEFFICIENTS
C  --------------------------------
 
      IF(J.EQ.2) THEN
         DO K=1,NCR
         DO I=1,NCF
         QLN(INX(I)) = COF(I,K)
         ENDDO
         DO I=1,NCF
         COF(I,K) = QLN(I)
         ENDDO
         ENDDO
      ENDIF
 
C  MAKE QLN FOR A SPECIFIC LATITUDE
C  --------------------------------
 
      SINLAT = COS(COLR(J))
 
      DO N=1,JCAP1
      N1 = 2*N-1
      N2 = N1+JCAP1X2
      QLN(N1  ) = PLN(N,1)
      QLN(N1+1) = PLN(N,1)
      QLN(N2  ) = PLN(N,2)
      QLN(N2+1) = PLN(N,2)
      ENDDO
 
      LP0 = 2*JCAP1X2-2*N0
      LP1 = JCAP1X2
      LP2 = 0
 
      DO N=1,JCAP
      LEN = JCAP1X2-2*(N+N0)
      DO L=1,LEN
      QLN(LP0+L) = (SINLAT*QLN(LP1+L)-DEP(LP1+L)*QLN(LP2+L))/DEP(LP0+L)
      ENDDO
      LP2 = LP1
      LP1 = LP0
      LP0 = LP0 + LEN
      ENDDO
 
 
C  SUM THE COEFFICIENTS FOR THIS LATITUDE
C  --------------------------------------
 
      LL   = 0
      FST  = 0
      HEM  = 1
      LEN  = JCAP1X2
      AHEM = SIGN(1.,XLAT(J))
 
      DO N=N0,JCAP1
      DO K=1,NCR
      DO L=1,LEN
      GRD(L,K) = COF(L+LL,K)*QLN(L+LL)*HEM + GRD(L,K)
      ENDDO
      ENDDO
 
      LL  = LL+LEN
      LEN = JCAP1X2-2*N
      HEM = HEM*AHEM
      FST = 1
      ENDDO
 
C     IF(NCF.EQ.MDIMC) THEN
C     WRITE(50)(COF(I,1),I=1,NCF)
C     WRITE(50)(QLN(I),I=1,NCF)
C     WRITE(50)(GRD(I,1),I=1,JCAP1X2)
C     ENDIF
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POLES(GRD,IQ,SNL,CSL,RSC)
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
 
C-CRA DIMENSION GRD(IMAX,JMAX,2),SNL(IMAX),CSL(IMAX),RSC(JMAX)
 
      DIMENSION GRD( 256 , 129 ,2)
      DIMENSION SNL( 256 ),CSL( 256 ),RSC( 129 )
 
      DATA PI180/.0174532/
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      RIMAX = 1./IMAX
 
C  SURFACE PRESSURE AND WINDS GET SPECIAL TREATMENT
C  ------------------------------------------------
 
      IF(IQ.EQ.2) THEN
         DO J=2,JMAX-1
         DO I=1,IMAX
         GRD(I,J,1) = 10.*EXP(GRD(I,J,1))
         ENDDO
         ENDDO
      ELSEIF(IQ.EQ.4) THEN
         DO J=2,JMAX-1
         DO I=1,IMAX
         GRD(I,J,1) = GRD(I,J,1)*RSC(J)
         GRD(I,J,2) = GRD(I,J,2)*RSC(J)
         ENDDO
         ENDDO
         GOTO 20
      ENDIF
 
C  AVERAGE THE NEAREST ZONE FOR SCALAR POLES
C  -----------------------------------------
 
10    DO J=1,JMAX,JMAX-1
      IF(J.EQ.JMAX) JN = JMAX-1
      IF(J.EQ.1   ) JN = 2
      POLE = 0
      DO I=1,IMAX
      POLE = POLE + GRD(I,JN,1)*RIMAX
      ENDDO
      DO I=1,IMAX
      GRD(I,J,1) = POLE
      ENDDO
      ENDDO
 
      RETURN
 
C  AVERAGE THE NEAREST ZONE FOR VECTOR POLES
C  -----------------------------------------
 
20    DO J=1,JMAX,JMAX-1
      IF(J.EQ.JMAX) JN = JMAX-1
      IF(J.EQ.1   ) JN = 2
      UPOLE = 0
      VPOLE = 0
      DO I=1,IMAX
      SNL(I) = -SNL(I)
      UPOLE = UPOLE + (GRD(I,JN,1)*CSL(I)-GRD(I,JN,2)*SNL(I))*RIMAX
      VPOLE = VPOLE + (GRD(I,JN,1)*SNL(I)+GRD(I,JN,2)*CSL(I))*RIMAX
      ENDDO
      DO I=1,IMAX
      GRD(I,J,1) =  UPOLE*CSL(I) + VPOLE*SNL(I)
      GRD(I,J,2) = -UPOLE*SNL(I) + VPOLE*CSL(I)
      ENDDO
      ENDDO
 
      RETURN
      END
      SUBROUTINE MAXMIN(F,IDIM,JDIM,IMAX,JMAX,KMAX)
C
      DIMENSION F(IDIM,JDIM,KMAX)
C
      DO 10 K=1,KMAX
C
      FMAX=F(1,1,K)
      FMIN=F(1,1,K)
C
      DO 20 J=1,JMAX
      DO 20 I=1,IMAX
      IF(FMAX.LE.F(I,J,K)) THEN
      FMAX=F(I,J,K)
      IIMAX=I
      JJMAX=J
      ENDIF
      IF(FMIN.GE.F(I,J,K)) THEN
      FMIN=F(I,J,K)
      IIMIN=I
      JJMIN=J
      ENDIF
   20 CONTINUE
C
      WRITE(6,100) K,FMAX,IIMAX,JJMAX,FMIN,IIMIN,JJMIN
  100 FORMAT(2X,'LEVEL=',I2,' MAX=',E10.4,' AT I=',I5,' J=',I5,
     1                      ' MIN=',E10.4,' AT I=',I5,' J=',I5)
C
   10 CONTINUE
C
      RETURN
      END
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C MAIN PROGRAM:  PREVENTS    PREPROCESSING EVENTS PROGRAM
C   PRGMMR: WOOLLEN          ORG: W/NMC2      DATE: 94-09-06
C
C ABSTRACT: WRITES BACKGROUND EVENTS INTO THE PREPBUFR DATA FILE IN
C   PREPARATION FOR QUALITY CONTROL AND ANALYSIS PROCESSING. EVENTS
C   WRITTEN INCLUDE FORECAST FIRST GUESS VALUES INTERPOLATED FROM
C   THE SPECTRAL SIGMA GUESS FILE, AND OBSERVATION ERRORS FROM THE
C   SSI ERROR SPECIFICATION FILE. FOR CASES WHERE THE OB ERROR IS
C   MISSING IN THE SSI ERROR FILE, AN EVENT IS RECORDED WITH A
C   QUALITY MARK OF 9, WHICH SCREENS THE OB FROM FURTHER PROCESSING
C   IN THE QC AND ANALYSIS STEPS. ALL SURFACE PRESSURE OBS ARE ALSO
C   CHECKED AGAINST THE MODEL SURFACE PRESSURE, WITH THOSE REPORTING
C   DIFFERENCES GREATER THAN 100 MB MARKED WITH A QUALITY OF 8. FINALLY,
C   SOME GENERAL RULES ARE APPLIED WHICH ARBITRARILY SCREEN CERTAIN
C   CLASSES OF OBSERVATIONS FROM FURTHER PROCESSING (SEE SUBROUTINE
C   FILTAN).
C
C PROGRAM HISTORY LOG:
C   94-01-06  J. WOOLLEN  ORIGINAL VERSION FOR REANALYSIS
C   94-09-06  J. WOOLLEN  VERSION FOR IMPLEMENTATION IN GBL SYSTEM
C   12-06-08  W. EBISUZAI added OPEN(UNIT=LUGES,
C                     form='unformatted), for ifort compiler
C
C USAGE:
C   INPUT FILES:
C     FORT.11  - PREPDA     - OBSERVATION FILE
C     FORT.12  - SGES       - SPECTRAL SIGMA GUESS FILE
C     FORT.13  - SSIERR     - OBSERVATION ERROR FILE
C     FORT.14  - NMCDATE    - NMC DATE FILE
C
C   OUTPUT FILES:
C     FORT.50  - PREPDV      - PROCESSED OBSERVATION FILE
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:    DATCHK,ETABLE,SETTERP,GETFC,GETOE,FILTAN,HTERP
C
C   LIBRARY:
C     COMMON:    W3LIB,GESRES,BUFRLIB
C
C   EXIT STATES:
C     COND =   0 - SUCCESSFUL RUN
C
C     IF ANY DISASTORIUS PROBLEM SITUATIONS ARE DETECTED THE PROGRAM
C     WILL ABORT (RC=132) WITH A MESSAGE DETAILING THE CALAMITY.
C
C
C REMARKS:
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C   MACHINE:  CRAY
C
C$$$
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      PROGRAM PREVENTS
 
C-CRA TASKCOMMON /REPORT/ SID,XOB,YOB,DHR,TYP,NLEV,OBS(8,255),BAK(8,255)
      COMMON /REPORTI/ NLEV
      COMMON /REPORT/ SID,XOB,YOB,DHR,TYP,OBS(8,255),BAK(8,255)
 
      COMMON /PCODE / PVCD,VTCD
 
      CHARACTER*80 HEADR,OBSTR,FCSTR,OESTR
      CHARACTER*8  SUBSET
      DIMENSION    HDR(10),ZOQ(4)
 
      DATA HEADR /'SID XOB YOB DHR TYP ELV                   '/
      DATA OBSTR /'POB QOB TOB ZOB UOB VOB PWO CAT           '/
      DATA FCSTR /'PFC QFC TFC ZFC UFC VFC PWF               '/
      DATA OESTR /'POE QOE TOE ZOE WOE PWE                   '/
 
      DATA BMISS /10E10/
      DATA LUBFU /50/
      DATA LUBFA /51/
      DATA LUBFR /52/
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
 
C     CALL W3LOG('$S','$M')
 
      LUBFI = 11
      LUGES = 12
      LUERR = 13
      LUDAT = 14
 
      PRINT*
      PRINT*,'******  BEGINNING PREVENTS PROCESSING ******'
      PRINT*
 
C  OPEN I/O DATA FILES AND CHECK DATES
C  -----------------------------------
 
      CALL DATEBF(LUBFI,IYY,IMM,IDD,IHH,IDATE)
      CALL DATCHK(LUGES,LUDAT,IDATE)
 
C  OBTAIN THE FIRST GUESS
C  ----------------------
 
      CALL GESRES(LUGES,IGES)
      IF(IGES.NE.0) CALL SABORT('PREVENTS - UNABLE TO TRANSFORM GUESS')
 
C  READ ERROR FILES AND SET INTERPOLATION FACTORS
C  ----------------------------------------------
 
      CALL ETABLE(LUERR)
      CALL SETTERP
 
C  OPEN THE FILES AND GET THE QC CODES
C  -----------------------------------
 
      CALL OPENBF(LUBFI,'IN ',LUBFI)
      CALL OPENBF(LUBFU,'OUT',LUBFI)
      CALL OPENBF(LUBFA,'OUT',LUBFI)
      CALL OPENBF(LUBFR,'OUT',LUBFI)
 
      CALL UFBQCD(LUBFI,'PREVENT',PVCD)
      CALL UFBQCD(LUBFI,'VIRTMP ',VTCD)
 
C**---------------------------------------------------------------------
C**---------------------------------------------------------------------
 
C  LOOP THROUGH THE INPUT MESSAGES
C  -------------------------------
 
c      icnt=0
      DO WHILE(IREADMG(LUBFI,SUBSET,IDATE).EQ.0)
c        icnt=icnt+1
c	write(*,*) icnt,': subset=',subset
 
      IF(SUBSET.EQ.'ADPUPA') THEN
         LUBFO = LUBFU
      ELSEIF(SUBSET.EQ.'AIRCFT') THEN
         LUBFO = LUBFA
      ELSEIF(SUBSET.EQ.'AIRCAR') THEN
         LUBFO = LUBFA
      ELSE
         LUBFO = LUBFR
      ENDIF
      CALL OPENMB(LUBFO,SUBSET,IDATE)
      DO WHILE(IREADSB(LUBFI).EQ.0)
 
C  READ A REPORT
C  -------------
 
      CALL UFBINT(LUBFI,HDR,10,  1,NLEV,HEADR)
      CALL UFBINT(LUBFI,OBS, 8,255,NLEV,OBSTR)
      SID = HDR(1)
      XOB = HDR(2)
      YOB = HDR(3)
      DHR = HDR(4)
      TYP = HDR(5)
      ELV = HDR(6)
      CALL UFBCPY(LUBFI,LUBFO)
 
C  FILL IN MISSING SURFACE ELEVATION IF POSSIBLE
C  ---------------------------------------------
 
      CALL UFBINT(LUBFI,ZSF,1,1,NSF,'CAT=0 ZOB')
      IF(NSF.EQ.1 .AND. ELV.GE.BMISS .AND. ZSF.LT.BMISS) THEN
         CALL UFBINT(LUBFO,ZSF,1,1,NSF,'ELV')
      ELSEIF(NSF.EQ.1 .AND. ELV.LT.BMISS .AND. ZSF.GE.BMISS) THEN
         ZOQ(1) = ELV
         ZOQ(2) = 2
         ZOQ(3) = PVCD
         ZOQ(4) = 0
         CALL UFBINT(LUBFO,ZOQ,4,1,IRET,'CAT=0 ZOB ZQM ZPC ZRC')
         CALL UFBINT(-LUBFO,OBS,8,255,NLEV,OBSTR)
      ENDIF
 
C  MAKE THE PREVENTS AND VIRTUAL TEMPERATURE EVENTS
C  ------------------------------------------------
 
      CALL GETFC
      CALL UFBINT(LUBFO,BAK,8,NLEV,IRET,FCSTR)
      CALL GETOE
      CALL UFBINT(LUBFO,BAK,8,NLEV,IRET,OESTR)
      CALL FILTAN(LUBFO)
 
      IF(SUBSET.NE.'ADPUPA' .AND. SUBSET.NE.'SATEMP') THEN
         CALL VTPEVN(LUBFI,LUBFO)
      ENDIF
 
C  END OF READ LOOPS - WRITE THE PREVENTED SUBSET
C  ----------------------------------------------
      CALL WRITSB(LUBFO)
      ENDDO
      ENDDO
      CALL CLOSMG(LUBFU)
      CALL CLOSMG(LUBFA)
      CALL CLOSMG(LUBFR)
 
C**---------------------------------------------------------------------
C**---------------------------------------------------------------------
 
C  CLOSE THE BUFR FILES
C  --------------------
 
      CALL CLOSBF(LUBFI)
      CALL CLOSBF(LUBFU)
      CALL CLOSBF(LUBFA)
      CALL CLOSBF(LUBFR)
 
C  END OF PREVENTION
C  -----------------
 
C     CALL W3LOG('$E')
 
      STOP
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DATCHK(LUGES,LUDAT,IDATE)
 
      CHARACTER*8  O85LAB(4),DATE
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  SIGMA FILE
C  ----------
 
C wne 6-2012      REWIND LUGES
      open(unit=LUGES,form='unformatted')
      READ(LUGES,END=900,ERR=900) O85LAB
C     READ(LUGES,END=900,ERR=900) FHOUR
      READ(LUGES,END=900,ERR=900) FHOUR,IH,IM,ID,IY
	write(*,*) 'prevents: sigges - fhour,idate(4)=', FHOUR,IH,IM,ID,IY
      iy = mod(iy,100)
      JFHOUR = FHOUR
C     CALL W3FS11(O85LAB(2),IY,IM,ID,IH,1)
      if (iy.ne.0) then
         CALL W3FS03(IDI,IH,IY,IM,ID,0)
         CALL W3FS15(IDI,JFHOUR,IDV)
         CALL W3FS03(IDV,IH,IY,IM,ID,1)
      else
         iy=iy+80
         CALL W3FS03(IDI,IH,IY,IM,ID,0)
         CALL W3FS15(IDI,JFHOUR,IDV)
         CALL W3FS03(IDV,IH,IY,IM,ID,1)
         iy=iy-80
         iy=mod(iy,100)
      endif
      WRITE(DATE,'(4I2)') IY,IM,ID,IH
      DO I=1,8
      IF(DATE(I:I).EQ.' ') DATE(I:I) = '0'
      ENDDO
      PRINT'(''GUESS VALID AT  '',A8)',DATE
 
C  PREPDA FILE
C  -----------
 
      WRITE(DATE,'( I8)') IDATE
      READ( DATE,'(4I2)') JY,JM,JD,JH
      DO I=1,8
      IF(DATE(I:I).EQ.' ') DATE(I:I) = '0'
      ENDDO
      PRINT'(''DATA  VALID AT  '',A8)',DATE
 
C  NMCDATE FILE
C  ------------
 
      REWIND LUDAT
      READ(LUDAT,'(8X,4I2)',END=901,ERR=901) KY,KM,KD,KH
      WRITE(DATE,'(4I2)') KY,KM,KD,KH
      DO I=1,8
      IF(DATE(I:I).EQ.' ') DATE(I:I) = '0'
      ENDDO
      PRINT'(''NMCDATE         '',A8)',DATE
 
C  VALID DATES MUST MATCH
C  ----------------------
 
      IF(IY.NE.JY .OR. IM.NE.JM .OR. ID.NE.JD .OR. IH.NE.JH .OR.
     .   IY.NE.KY .OR. IM.NE.KM .OR. ID.NE.KD .OR. IH.NE.KH)
     .   GOTO 902
 
      RETURN
900   CALL SABORT('DATCHK - BAD OR MISSING GUESS FILE     ')
901   CALL SABORT('DATCHK - BAD OR MISSING NMCDATE FILE   ')
902   CALL SABORT('DATCHK - DATES DONT MATCH              ')
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ETABLE(LUNIT)
 
      COMMON /ETRP/ ERRS(300,33,6)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  READ THE OBSERVATION ERROR TABLES
C  ---------------------------------
 
      REWIND LUNIT
C-CRA ERRS = 0
C     COMMON /ETRP/ ERRS(300,33,6)
      DO I=1,300*33*6
        ERRS(I,1,1) = 0
      ENDDO
 
10    READ(LUNIT,'(1X,I3)',END=100) KX
      DO K=1,33
      READ(LUNIT,'(1X,6E12.5)') (ERRS(KX,K,M),M=1,6)
      ENDDO
      GOTO 10
 
100   RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FILTAN(LUBFO)
 
C-CRA TASKCOMMON /REPORT/ SID,XOB,YOB,DHR,TYP,NLEV,OBS(8,255),BAK(8,255)
C-CRA TASKCOMMON /GUESS / PS,ZS,T(100),U(100),V(100),Q(100)
      COMMON /REPORTI/ NLEV
      COMMON /REPORT/ SID,XOB,YOB,DHR,TYP,OBS(8,255),BAK(8,255)
      COMMON /GUESS / PS,ZS,T(100),U(100),V(100),Q(100)
      COMMON /PCODE / PVCD,VTCD
 
      CHARACTER*40 PEVN,QEVN,TEVN,WEVN,PWVN
      DIMENSION    PEV(4,255),QEV(4,255),TEV(4,255),WEV(5,255)
      DIMENSION    PWV(4,255)
      LOGICAL      ACARS,SATEMP,SOLN60,SOLS60,IOLN60,IOLS60,SSMI,MISS
      LOGICAL      DN2FAR,REJP,REJPS,REJT,REJQ,REJW,REJPW
 
      DATA PEVN /'POB PQM PPC PRC     '/
      DATA QEVN /'QOB QQM QPC QRC     '/
      DATA TEVN /'TOB TQM TPC TRC     '/
      DATA WEVN /'UOB VOB WQM WPC WRC '/
      DATA PWVN /'PWO PWQ PWP PWR     '/
 
      DATA BMISS /10E10/
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  LOGICAL SWITCHES FOR OBSERVATION LOCATION FILTERING
C  ---------------------------------------------------
 
      SSMI   = TYP.EQ.283
      ACARS  = TYP.EQ.233
      SATEMP = TYP.GE.160.AND.TYP.LE.179
      SOLN60 = TYP.GE.160.AND.TYP.LE.163.AND.YOB.GE.-60
      IOLN60 = TYP.GE.165.AND.TYP.LE.168.AND.YOB.GE.-60
      SOLS60 = (TYP.EQ.160.OR.TYP.EQ.162.OR.TYP.EQ.163).AND.YOB.LT.-60
      IOLS60 = (TYP.EQ.165.OR.TYP.EQ.167.OR.TYP.EQ.168).AND.YOB.LT.-60
 
C  CLEAR THE EVENT ARRAYS
C  ----------------------
 
C-CRA PEV = BMISS
C-CRA QEV = BMISS
C-CRA TEV = BMISS
C-CRA WEV = BMISS
C-CRA PWV = BMISS
C     DIMENSION    PEV(4,255),QEV(4,255),TEV(4,255),WEV(5,255)
      DO I=1,4*255
        PEV(I,1) = BMISS
        QEV(I,1) = BMISS
        TEV(I,1)= BMISS
      ENDDO
      DO I=1,5*255
      WEV(I,1) = BMISS
      ENDDO
      DO I=1,4*255
        PWV(I,1) = BMISS
      ENDDO
 
      MAXPEV = 0
      MAXQEV = 0
      MAXTEV = 0
      MAXWEV = 0
      MAXPWV = 0
 
C  LOOP OVER LEVEL APPLYING UNDERGROUND FILTERING AND SPECIAL RULES
C  ----------------------------------------------------------------
 
      DO 50 L=1,NLEV
 
      POB = OBS(1,L)
      QOB = OBS(2,L)
      TOB = OBS(3,L)
      ZOB = OBS(4,L)
      UOB = OBS(5,L)
      VOB = OBS(6,L)
      PWO = OBS(7,L)
      CAT = OBS(8,L)
 
      DN2FAR = NLEV.EQ.1 .AND. TYP.EQ.120
      REJ = 9
 
C  ANY PRESSURE MORE THAN 100 MB BELOW MODEL SURFACE IS REJECTED
C  -------------------------------------------------------------
C  ANY ZERO OR NEGATIVE PRESSURE IS REJECTED
C  -----------------------------------------
 
      IF(POB.LT.BMISS .AND. (POB-PS.GE.100. .OR. POB.LE.0.)) THEN
         REJ = 8
         DN2FAR = .TRUE.
         PEV(1,L) = POB
         PEV(2,L) = REJ
         PEV(3,L) = PVCD
         PEV(4,L) = REJ
         MAXPEV = L
      ENDIF
 
C  RULES FOR SURFACE PRESSURE
C  --------------------------
 
      IF(POB.LT.BMISS .AND. CAT.EQ.0) THEN
         REJPS = OEF(POB,TYP,5).GE.BMISS         .OR.
     .           ABS(POB-PS).GE.100.             .OR.
     .           POB.LE.450.                     .OR.
     .           POB.GE.1100.
         IF(REJPS.OR.DN2FAR) THEN
            DN2FAR = .TRUE.
            PEV(1,L) = POB
            PEV(2,L) = REJ
            PEV(3,L) = PVCD
            PEV(4,L) = REJ
            MAXPEV = L
         ENDIF
      ENDIF
 
C  RULES FOR TEMPERATURE
C  ---------------------
 
      IF(TOB.LT.BMISS) THEN
         REJT = OEF(POB,TYP,2).GE.BMISS          .OR.
     .          (SOLN60 .AND. NINT(POB).GE.100)  .OR.
     .          (IOLN60 .AND. NINT(POB).GE.100)  .OR.
     .          (SOLS60 .AND. NINT(POB).GT.100)  .OR.
     .          (IOLS60 .AND. NINT(POB).GT.100)
         IF(REJT.OR.DN2FAR) THEN
            TEV(1,L) = TOB
            TEV(2,L) = REJ
            TEV(3,L) = PVCD
            TEV(4,L) = REJ
            MAXTEV = L
         ENDIF
      ENDIF
 
C  RULES FOR SPECIFIC HUMIDITY
C  ---------------------------
 
      IF(QOB.LT.BMISS) THEN
         REJQ = OEF(POB,TYP,3).GE.BMISS          .OR.
     .          TOB.GE.BMISS                     .OR.
     .          TOB.LE.-150.                     .OR.
     .          POB.LE.300                       .OR.
     .          QOB.LE.0.                        .OR.
     .          SATEMP                           .OR.
     .          REJT
         IF(REJQ.OR.DN2FAR) THEN
            QEV(1,L) = QOB
            QEV(2,L) = REJ
            QEV(3,L) = PVCD
            QEV(4,L) = REJ
            MAXQEV = L
         ENDIF
      ENDIF
 
C  RULES FOR WINDS
C  ---------------
 
      IF(UOB.LT.BMISS .OR. VOB.LT.BMISS) THEN
         REJW = OEF(POB,TYP,4).GE.BMISS .OR. (ACARS .AND. POB.GT.700)
         MISS = (UOB.GE.BMISS.OR.VOB.GE.BMISS) .AND. .NOT.SSMI
         IF(REJW.OR.MISS.OR.DN2FAR) THEN
            WEV(1,L) = UOB
            WEV(2,L) = VOB
            WEV(3,L) = REJ
            WEV(4,L) = PVCD
            WEV(5,L) = REJ
            MAXWEV = L
         ENDIF
      ENDIF
 
C  RULES FOR PRECIPITABLE WATER
C  ----------------------------
 
      IF(PWO.LT.BMISS) THEN
         REJPW = OEF(POB,TYP,6).GE.BMISS    .OR.
     .           POB.LE.20
         IF(REJPW) THEN
            PWV(1,L) = PWO
            PWV(2,L) = 9
            PWV(3,L) = PVCD
            PWV(4,L) = 9
            MAXPWV = L
         ENDIF
      ENDIF
 
50    ENDDO
 
C  APPLY THE PROPER EVENTS
C  -----------------------
 
      IF(MAXPEV.GT.0) CALL UFBINT(LUBFO,PEV,4,MAXPEV,IRET,PEVN)
      IF(MAXQEV.GT.0) CALL UFBINT(LUBFO,QEV,4,MAXQEV,IRET,QEVN)
      IF(MAXTEV.GT.0) CALL UFBINT(LUBFO,TEV,4,MAXTEV,IRET,TEVN)
      IF(MAXWEV.GT.0) CALL UFBINT(LUBFO,WEV,5,MAXWEV,IRET,WEVN)
      IF(MAXPWV.GT.0) CALL UFBINT(LUBFO,PWV,4,MAXPWV,IRET,PWVN)
 
      RETURN
      END
C-----------------------------------------------------------------------
C  GETFC - INTERPOLATE FIRST GUESS TO OB LOCATIONS
C-----------------------------------------------------------------------
      SUBROUTINE GETFC
 
C-CRA TASKCOMMON /REPORT/ SID,XOB,YOB,DHR,TYP,NLEV,OBS(8,255),BAK(8,255)
      COMMON /REPORTI/ NLEV
      COMMON /REPORT/ SID,XOB,YOB,DHR,TYP,OBS(8,255),BAK(8,255)
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
      COMMON /SIGFACI/ MAXPRS
      COMMON /SIGFAC/ SITERP(1200),SLTERP(1200)
C-CRA TASKCOMMON /GUESS / PS,ZS,T(100),U(100),V(100),Q(100)
      COMMON /GUESS / PS,ZS,T(100),U(100),V(100),Q(100)
 
      PARAMETER(KDIM=500)
      DIMENSION PINT(KDIM),ZINT(KDIM)
 
      DATA BMISS / 10E10  /
      DATA TZERO / 273.15 /
      DATA BETAP / .0552  /
      DATA BETA  / .00650 /
      DATA ROG   / 29.261 /
      DATA G     / 9.81   /
      DATA R     / 287.05 /
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      IF(KMAX.GT.KDIM) THEN
        PRINT *,'KMAX.GT.KDIM    KMAX=',KMAX
        CALL SABORT
      ENDIF
 
C  CLEAR THE BACKGROUND EVENT ARRAY
C  --------------------------------
 
      DO I=1,8*255
        BAK(I,1) = BMISS
      ENDDO
 
C  GET SIGMA PROFILE AT OB LOCATION
C  --------------------------------
 
      CALL HTERP(XOB,YOB)
      PS000 = 1000./PS
      PSIG  = PS*SL(1)
      PINT(1) = PS
      ZINT(1) = ZS
      DO K=2,KMAX
      K0 = K-1
      PINT(K) = PS*SI(K)
      ZINT(K) = ZINT(K0) - ROG*T(K0)*LOG(PINT(K)/PINT(K0))
      ENDDO
 
C  INTERPOLATE GUESS PROFILES TO OB PRESSURES
C  ------------------------------------------
 
      DO 10 L=1,NLEV
 
      POB = OBS(1,L)
      QOB = OBS(2,L)
      TOB = OBS(3,L)
      ZOB = OBS(4,L)
      UOB = OBS(5,L)
      VOB = OBS(6,L)
      PWO = OBS(7,L)
      CAT = OBS(8,L)
      IF(POB.LE.0. .OR. POB.GE.BMISS) GOTO 10
 
      IP  = POB*PS000
      IP  = MIN(IP,MAXPRS)
      IP  = MAX(IP,1)
 
C  SURFACE PRESSURE
C  ----------------
 
      IF(CAT.EQ.0 .AND. ZOB.LT.BMISS) THEN
         TS = T(1) + (PS-PSIG)*BETAP
         DZ  = ZOB-ZS
         TM  = TS - DZ*BETA*.5
         PFC = PS*EXP(-DZ/(TM*ROG))
      ELSE
         PFC = BMISS
      ENDIF
 
C  SPECIFIC HUMIDITY
C  -----------------
 
      IF(QOB.LT.BMISS) THEN
         LB = SLTERP(IP)
         WT = SLTERP(IP)-LB
         LA = MIN(LB+1,KMAX)
         QOB = Q(LB) + (Q(LA)-Q(LB))*WT
      ENDIF
 
C  TEMPERATURE
C  -----------
 
      IF(TOB.LT.BMISS) THEN
         IF(POB.GT.PSIG) THEN
            TOB = T(1) + (POB-PSIG)*BETAP
         ELSE
            LB = SLTERP(IP)
            WT = SLTERP(IP)-LB
            LA = MIN(LB+1,KMAX)
            TOB = T(LB) + (T(LA)-T(LB))*WT
         ENDIF
         TOB = TOB - TZERO
      ENDIF
 
C  HEIGHT
C  ------
 
      IF(ZOB.LT.BMISS) THEN
         IF(POB.GT.PSIG) THEN
            TM = T(1) + (.5*(PINT(1)+POB)-PSIG)*BETAP
            ZOB = ZINT(1) - ROG*TM*LOG(POB/PINT(1))
         ELSE
            LI = SITERP(IP)
            MP = (POB+PINT(LI))*PS000*.5
            MP = MAX(MIN(MP,MAXPRS),1)
            LB = SLTERP(MP)
            WT = SLTERP(MP)-LB
            LA = MAX(MIN(LB+1,KMAX),1)
            TM = T(LB) + (T(LA)-T(LB))*WT
            ZOB = ZINT(LI) - ROG*TM*LOG(POB/PINT(LI))
         ENDIF
      ENDIF
 
C  U AND V COMPONENTS
C  ------------------
 
      IF(UOB.LT.BMISS .OR. VOB.LT.BMISS) THEN
         LB = SLTERP(IP)
         WT = SLTERP(IP)-LB
         LA = MIN(LB+1,KMAX)
         UOB = U(LB) + (U(LA)-U(LB))*WT
         VOB = V(LB) + (V(LA)-V(LB))*WT
      ENDIF
 
C  PRECIPITABLE WATER
C  ------------------
 
      PWO = BMISS
 
C  RELATIVE HUMIDITY
C  -----------------
 
      RHO = BMISS
 
C  SCATTER THE PROPER FORECAST VALUES
C  ----------------------------------
 
      BAK(1,L) = PFC
      BAK(2,L) = QOB
      BAK(3,L) = TOB
      BAK(4,L) = ZOB
      BAK(5,L) = UOB
      BAK(6,L) = VOB
      BAK(7,L) = PWO
      BAK(8,L) = RHO
 
10    ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE GETOE
 
C-CRA TASKCOMMON /REPORT/ SID,XOB,YOB,DHR,TYP,NLEV,OBS(8,255),BAK(8,255)
      COMMON /REPORTI/ NLEV
      COMMON /REPORT/ SID,XOB,YOB,DHR,TYP,OBS(8,255),BAK(8,255)
 
      DATA BMISS /10E10/
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  CLEAR THE EVENT ARRAY
C  ---------------------
 
      DO I=1,8*255
        BAK(I,1) = BMISS
      ENDDO
 
C  LOOP OVER LEVELS LOOKING UP THE OBSERVATION ERROR
C  -------------------------------------------------
 
      DO L=1,NLEV
 
      POB = OBS(1,L)
      QOB = OBS(2,L)
      TOB = OBS(3,L)
      ZOB = OBS(4,L)
      WOB = MAX(OBS(5,L),OBS(6,L))
      PWO = OBS(7,L)
      CAT = OBS(8,L)
 
      IF(CAT.EQ.0    ) BAK(1,L) = OEF(POB,TYP,5)
      IF(QOB.LT.BMISS) BAK(2,L) = OEF(POB,TYP,3)
      IF(TOB.LT.BMISS) BAK(3,L) = OEF(POB,TYP,2)
      IF(WOB.LT.BMISS) BAK(5,L) = OEF(POB,TYP,4)
      IF(PWO.LT.BMISS) BAK(6,L) = OEF(POB,TYP,6)
 
      ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE GUSER - GESRES USER INTERFACE FOR PREPFIT (PS,ZS,T,U,V)
C-----------------------------------------------------------------------
      SUBROUTINE GUSER(GRD,IQ,LEV)
 
      COMMON /SIGMASI/IMAX,JMAX,KMAX
      COMMON /SIGMAS/DLAT,DLON,SL(100),SI(101)
 
C-CRA DIMENSION   GRD(IMAX,JMAX,2)
      DIMENSION   GRD( 256 , 129 ,2)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  PACK 2D GUESS FIELD INTO BIT PACKED ARRAYS
C  ------------------------------------------
 
      IF(IQ.EQ.1) THEN
         DO J=1,JMAX
         DO I=1,IMAX
         CALL ZSP(I,J,GRD(I,J,1))
         ENDDO
         ENDDO
      ELSEIF(IQ.EQ.2) THEN
         DO J=1,JMAX
         DO I=1,IMAX
         CALL PSP(I,J,GRD(I,J,1))
         ENDDO
         ENDDO
      ELSEIF(IQ.EQ.3) THEN
         DO J=1,JMAX
         DO I=1,IMAX
         CALL TP(I,J,LEV,GRD(I,J,1))
         ENDDO
         ENDDO
      ELSEIF(IQ.EQ.4) THEN
         DO J=1,JMAX
         DO I=1,IMAX
         CALL UP(I,J,LEV,GRD(I,J,1))
         CALL VP(I,J,LEV,GRD(I,J,2))
         ENDDO
         ENDDO
      ELSEIF(IQ.EQ.5) THEN
         DO J=1,JMAX
         DO I=1,IMAX
         CALL QP(I,J,LEV,MAX(0.,GRD(I,J,1))*1E6)
         ENDDO
         ENDDO
      ENDIF
 
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE HTERP - 2D LINEAR HORIZONTAL INTERPOLATION
C-----------------------------------------------------------------------
      SUBROUTINE HTERP(XOB,YOB)
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
C-CRA TASKCOMMON /GUESS / PSI,ZSI,TI(100),UI(100),VI(100),QI(100)
      COMMON /GUESS / PSI,ZSI,TI(100),UI(100),VI(100),QI(100)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
C  CALCULATE HORIZONTAL WEIGHTS AND INTERPOLATE
C  --------------------------------------------
 
      WX = XOB/DLON + 1.0
      I0 = WX
      I1 = MOD(I0,IMAX) + 1
      WX = WX-I0
 
      WY = (YOB+90.)/DLAT + 1.0
      J0 = WY
      J1 = MIN(J0+1,JMAX)
      WY = WY-J0
 
C  HTERP FOR SURFACE HEIGHT
C  ------------------------
 
      P1  = ZS(I0,J0)
      P2  = ZS(I0,J1)
      P3  = ZS(I1,J0)
      P4  = ZS(I1,J1)
      P5  = P1+(P2-P1)*WY
      P6  = P3+(P4-P3)*WY
      ZSI = P5+(P6-P5)*WX
 
C  HTERP FOR SURFACE PRESSURE
C  --------------------------
 
      P1  = PS(I0,J0)
      P2  = PS(I0,J1)
      P3  = PS(I1,J0)
      P4  = PS(I1,J1)
      P5  = P1+(P2-P1)*WY
      P6  = P3+(P4-P3)*WY
      PSI = P5+(P6-P5)*WX
 
C  HTERP FOR UPA T,U,V,Q
C  ---------------------
 
      DO K=1,KMAX
 
      P1 = T(I0,J0,K)
      P2 = T(I0,J1,K)
      P3 = T(I1,J0,K)
      P4 = T(I1,J1,K)
      P5 = P1+(P2-P1)*WY
      P6 = P3+(P4-P3)*WY
      TI(K) = P5+(P6-P5)*WX
 
      P1 = U(I0,J0,K)
      P2 = U(I0,J1,K)
      P3 = U(I1,J0,K)
      P4 = U(I1,J1,K)
      P5 = P1+(P2-P1)*WY
      P6 = P3+(P4-P3)*WY
      UI(K) = P5+(P6-P5)*WX
 
      P1 = V(I0,J0,K)
      P2 = V(I0,J1,K)
      P3 = V(I1,J0,K)
      P4 = V(I1,J1,K)
      P5 = P1+(P2-P1)*WY
      P6 = P3+(P4-P3)*WY
      VI(K) = P5+(P6-P5)*WX
 
      P1 = Q(I0,J0,K)
      P2 = Q(I0,J1,K)
      P3 = Q(I1,J0,K)
      P4 = Q(I1,J1,K)
      P5 = P1+(P2-P1)*WY
      P6 = P3+(P4-P3)*WY
      QI(K) = P5+(P6-P5)*WX
 
      ENDDO
 
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION OEF(P,TYP,IE)
 
      COMMON /ETRP/ERRS(300,33,6)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      OEF = 10E10
      KX  = TYP
 
C  LOOK UP ERRORS FOR PARTICULAR OB TYPES
C  --------------------------------------
 
      IF(IE.GE.2 .AND. IE.LE.4) THEN
         DO LA=1,33
         IF(P.GE.ERRS(KX,LA,1)) GOTO 10
         ENDDO
10       LB = LA-1
         IF(LB.EQ.33) LA = 6
         IF(LB.EQ.33) LB = 5
         IF(LB.EQ. 0) THEN
            OEF = ERRS(KX,1,IE)
         ELSE
            DEL = (P-ERRS(KX,LB,1))/(ERRS(KX,LA,1)-ERRS(KX,LB,1))
            OEF = (1.-DEL)*ERRS(KX,LB,IE) + DEL*ERRS(KX,LA,IE)
         ENDIF
      ELSEIF(IE.EQ.5) THEN
         OEF = ERRS(KX,1,5)
      ELSEIF(IE.EQ.6) THEN
         OEF = ERRS(KX,1,6)
      ENDIF
 
C  SET MISSING ERROR VALUE TO 10E10
C  --------------------------------
 
      IF(OEF.GE.5E5) OEF = 10E10
 
      RETURN
      END
C-----------------------------------------------------------------------
C  SUBROUTINE SETTERP
C-----------------------------------------------------------------------
      SUBROUTINE SETTERP
 
      COMMON /SIGMASI/ IMAX,JMAX,KMAX
      COMMON /SIGMAS/ DLAT,DLON,SL(100),SI(101)
      COMMON /SIGFACI/ MAXPRS
      COMMON /SIGFAC/ SITERP(1200),SLTERP(1200)
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      MAXPRS  = 1200
 
C  COMPUTE TABLE FOR THE SIGMA MIDPOINTS RELATIVE TO 1000 MB
C  ---------------------------------------------------------
 
      DO 20 I=1,MAXPRS
      P = I
      DO 10 LA=1,KMAX
      IF(P.GE.SL(LA)*1000.) GOTO 15
10    CONTINUE
15    IF(LA.GT.KMAX) THEN
         LA = KMAX
         LB = KMAX
      ELSE IF(LA.EQ.1) THEN
         LA = 1
         LB = 1
      ELSE
         LB = LA-1
      ENDIF
      IF(LA.EQ.LB) THEN
         WP = 0.
      ELSE
         PA = SL(LA)*1000.
         PB = SL(LB)*1000.
         WP = LOG(P/PB)/LOG(PA/PB)
      ENDIF
      SLTERP(I) = LB + WP
20    CONTINUE
 
C  COMPUTE TABLE FOR THE SIGMA INTERFACES RELATIVE TO 1000 MB
C  ----------------------------------------------------------
 
      DO 40 I=1,MAXPRS
      P = I
      DO 30 LA=1,KMAX
      IF(P.GE.SI(LA)*1000.) GOTO 35
30    CONTINUE
35    SITERP(I) = MAX(LA-1,1)
40    CONTINUE
 
      RETURN
      END
C-----------------------------------------------------------------------
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    VTPEVN      CALCULATE VIRTUAL TEMPERATURE
C   PRGMMR: J. Woollen       ORG: W/NMC20    DATE: 94-MM-DD
C
C ABSTRACT: CREATE VIRTUAL TEMPERATURE EVENTS WITHIN PREVENTS CODE FOR
C   ALL MESSAGE TYPES EXCEPT ADPUPA AND SATEMP
C
C PROGRAM HISTORY LOG:
C   95-05-17  J. WOOLLEN
C
C USAGE:    CALL VTPEVN(LUBFI,LUBFO,SUBSET)
C   INPUT ARGUMENT LIST:
C     LUBFI    - BUFR INPUT  FILE UNIT
C     LUBFO    - BUFR OUTPUT FILE UNIT
C     SUBSET   - MESSAGE TYPE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN77
C   MACHINE:  CRAY
C
C$$$
      SUBROUTINE VTPEVN(LUBFI,LUBFO)
 
C-CRA TASKCOMMON /REPORT/ SID,XOB,YOB,DHR,TYP,NLEV,OBS(8,255),BAK(8,255)
      COMMON /REPORTI/ NLEV
      COMMON /REPORT/ SID,XOB,YOB,DHR,TYP,OBS(8,255),BAK(8,255)
      COMMON /PCODE / PVCD,VTCD
 
      CHARACTER*80 EVNSTR
      DIMENSION    TDP(255),TQM(255),QQM(255)
      LOGICAL      EVN
 
      DATA EVNSTR /'QOB QQM QRC QPC TOB TQM TRC TPC'/
      DATA BMISS /10E10/
 
C-----------------------------------------------------------------------
C FCNS BELOW CONVERT TEMP/TD (K) & PRESS (MB) INTO SAT./ SPEC. HUM.(G/G)
C-----------------------------------------------------------------------
      ES(T) = 6.1078*EXP((17.269*(T - 273.16))/((T - 273.16)+237.3))
      QS(T,P) = (0.622*ES(T))/(P-(0.378*ES(T)))
C-----------------------------------------------------------------------
 
C  CLEAR TEMPERATURE AND Q EVENTS
C  ------------------------------
 
      EVN = .FALSE.
      DO I=1,8*255
        BAK(I,1) = BMISS
      ENDDO
 
C  GET DEWPOINT TEMPERATURE AND CURRENT T,Q QUALITY MARKS
C  ------------------------------------------------------
 
      CALL UFBINT(-LUBFO,TDP,1,255,NLTD,'TDO')
      CALL UFBINT(-LUBFO,TQM,1,255,NLTQ,'TQM')
      CALL UFBINT(-LUBFO,QQM,1,255,NLQQ,'QQM')
      IF(NLTD.NE.NLEV) CALL SABORT('VTPEVN - NLTD<>NLEV')
      IF(NLTQ.NE.NLEV) CALL SABORT('VTPEVN - NLTQ<>NLEV')
      IF(NLQQ.NE.NLEV) CALL SABORT('VTPEVN - NLQQ<>NLEV')
 
C  COMPUTE VIRTUAL TEMPERATURE AND SPECIFIC HUMIDITY USING REPORTED DEWP
C  ---------------------------------------------------------------------
 
      DO L=1,NLEV
      POB = OBS(1,L)
      TDO = TDP(L)
      TOB = OBS(3,L)
      IF(POB.LT.BMISS .AND. TOB.LT.BMISS
     .                .AND. TDO.LT.BMISS
     .                .AND. QQM(L).LT.4) THEN
         QOB = QS(TDO+273.16,POB)
         BAK(1,L) = QOB*1E6
         BAK(2,L) = QQM(L)
         BAK(3,L) = 0
         BAK(4,L) = VTCD
         BAK(5,L) = (TOB+273.16)*(1.+.61*QOB)-273.16
         BAK(6,L) = TQM(L)
         BAK(7,L) = 0
         BAK(8,L) = VTCD
         EVN = .TRUE.
      ENDIF
      ENDDO
 
C  WRITE THE EVENTS TO THE BUFR FILE
C  ---------------------------------
 
      IF(EVN) CALL UFBINT(LUBFO,BAK,8,NLEV,IRET,EVNSTR)
 
      RETURN
      END
       SUBROUTINE W3FS03(IDATE,IHOUR,IYEAR,MONTH,IDAY,NN)
C$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
C
C SUBPROGRAM: W3FS03         NMC DATE WORD PACKER AND UNPACKER
C   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-24
C
C ABSTRACT: OBTAINS THE COMPONENTS OF THE NMC DATE WORD (SEE NMC
C   O.N. 84 AND 85) OR GIVEN ITS COMPONENTS, FORMS AN NMC TYPE
C   DATE WORD. W3FS03 IS THE SAME AS W3FS11 EXCEPT FOR THE ORDER OF
C   THE PARAMETERS.
C
C PROGRAM HISTORY LOG:
C 87-03-24  R.E.JONES   CONVERT TO CYBER 205 FORTRAN 200
C 89-10-13  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
C
C USAGE:  CALL W3FS03 (IDATE,IHOUR,IYEAR,MONTH,IDAY,NN)
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IDATE  ARG LIST  LEFT 4 BYTES OF INTEGER WORD
C                      WORD. (7TH WORD OF DATA FIELD OR 3RD WORD OF
C                      THE ID TABLE OF A BINARY FILE IF IN HALF
C                      PRECISION ARRAY OR THE 4TH WORD OR 2ND WORD IF
C                      IN INTEGER ARRAY).
C     IHOUR  ARG LIST  HOUR
C     IYEAR  ARG LIST  YEAR (2 DIGITS)
C     MONTH  ARG LIST  MONTH
C     IDAY   ARG LIST  DAY
C     NN     ARG LIST  CODE:
C                       = 0  PACK IHOUR, IYEAR, MONTH, IDAY, INTO IDATE
C                      <> 0  UNPACK IDATE INTO IHOUR, IYEAR, MONTH, IDAY
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IDATE  ARG LIST  LEFT 4 BYTES OF INTEGER WORD
C                      WORD. (7TH WORD OF DATA FIELD OR 3RD WORD OF
C                      THE ID TABLE OF A BINARY FILE IF IN HALF
C                      PRECISION ARRAY OR THE 4TH WORD OR 2ND WORD IF
C                      IN INTEGER ARRAY).
C     IHOUR  ARG LIST  HOUR
C     IYEAR  ARG LIST  YEAR (2 DIGITS)
C     MONTH  ARG LIST  MONTH
C     IDAY   ARG LIST  DAY
C
C   SUBPRGRAMS CALLED:
C     NAMES   LIBRARY
C     ------------------------------------------------------- --------
C     CHAR mova2i                                              SYSTEM
C
C REMARKS:    WHEN NN.NE.0, THE INFORMATION IN IDATE MUST BE
C     FORMATTED AS DIAGRAMMED IN APPENDIX C OF NMC O.N. 84.
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CYAY Y-MP8/832
C
C$$$
C
      CHARACTER*1 IDATE(4)
C
      IF (NN.NE.0) THEN
C
        IYEAR = mova2i(IDATE(1))
        MONTH = mova2i(IDATE(2))
        IDAY  = mova2i(IDATE(3))
        IHOUR = mova2i(IDATE(4))
C
      ELSE
C
        IDATE(1) = CHAR(IYEAR)
        IDATE(2) = CHAR(MONTH)
        IDATE(3) = CHAR(IDAY)
        IDATE(4) = CHAR(IHOUR)
      ENDIF
C
      RETURN
      END
       SUBROUTINE W3FS11(IDATE,IYEAR,MONTH,IDAY,IHOUR,NN)
C$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
C
C SUBPROGRAM: W3FS11         NMC DATE WORD UNPACKER AND PACKER
C   AUTHOR: JONES,R.E.       ORG: W342       DATE: 87-03-24
C
C ABSTRACT: OBTAINS THE COMPONENTS OF THE NMC DATE WORD (NMC OFFICE
C   NOTES 84 AND 85), OR GIVEN ITS COMPONENTS, FORMS AN NMC TYPE DATE
C   WORD. W3FS11 IS THE SAME AS W3FS03 EXCEPT FOR THE ORDER OF THE
C   PARAMETERS.
C
C PROGRAM HISTORY LOG:
C   87-03-24  R.E.JONES   CONVERT TO CYBER 205 FORTRAN 200
C   89-10-13  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
C
C USAGE:  CALL W3FS11 (IDATE, IYEAR, MONTH, IDAY, IHOUR, NN)
C
C   INPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IDATE  ARG LIST  LEFT 4 BYTES OF INTEGER 64 BIT WORD, OR CAN BE
C                      CHARACTER*1 IDATE(4) OR CHARACTER*4 IDATE.
C                      IF OFFICE NOTE 85 LABEL USED AS 4 64 BIT WORDS,
C                      IDATE IS IN THE LEFT 32 BITS OF THE 2ND WORD.
C                      IF OFFICE NOTE 84 12 IDS. (6 64 BIT WORDS ON
C                      CRAY, DATE WORD IN LEFT 32 BITS OF 4TH CRAY ID
C                      WORD, OR 7TH 32 BIT ID WORD ON NAS.
C     IYEAR  ARG LIST  INTEGER   YEAR (2 DIGITS)
C     MONTH  ARG LIST  INTEGER   MONTH
C     IDAY   ARG LIST  INTEGER   DAY
C     IHOUR  ARG LIST  INTEGER   HOUR
C     NN     ARG LIST  INTEGER   CODE:
C                     .EQ. 0 PACK IYEAR, MONTH, IDAY, IHOUR INTO IDATE
C                     .NE. 0 UNPACK IDATE INTO IYEAR, MONTH, IDAY, IHOUR
C
C   OUTPUT VARIABLES:
C     NAMES  INTERFACE DESCRIPTION OF VARIABLES AND TYPES
C     ------ --------- -----------------------------------------------
C     IDATE  ARG LIST  LEFT 4 BYTES OF INTEGER 64 BIT WORD, OR CAN BE
C                      CHARACTER*1 IDATE(4) OR CHARACTER*4 IDATE.
C                      IF OFFICE NOTE 85 LABEL USED AS 4 64 BIT WORDS,
C                      IDATE IS IN THE LEFT 32 BITS OF THE 2ND WORD.
C                      IF OFFICE NOTE 84 12 IDS. (6 64 BIT WORDS ON
C                      CRAY, DATE WORD IN LEFT 32 BITS OF 4TH CRAY ID
C                      WORD, OR 7TH 32 BIT ID WORD ON NAS.
C     IYEAR  ARG LIST  INTEGER   YEAR (2 DIGITS)
C     MONTH  ARG LIST  INTEGER   MONTH
C     IDAY   ARG LIST  INTEGER   DAY
C     IHOUR  ARG LIST  INTEGER   HOUR
C
C   SUBROGRAMS CALLED:
C     NAMES                                                   LIBRARY
C     ------------------------------------------------------- --------
C     CHAR   mova2i                                            SYSTEM
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY Y-MP8/832
C
C$$$
C
      CHARACTER IDATE(4)
C
      IF (NN.NE.0) THEN
C
        IYEAR = mova2i(IDATE(1))
        MONTH = mova2i(IDATE(2))
        IDAY  = mova2i(IDATE(3))
        IHOUR = mova2i(IDATE(4))
C
      ELSE
C
        IDATE(1) = CHAR(IYEAR)
        IDATE(2) = CHAR(MONTH)
        IDATE(3) = CHAR(IDAY)
        IDATE(4) = CHAR(IHOUR)
      ENDIF
C
      RETURN
      END
      SUBROUTINE W3FS15(IDATE,JTAU,NDATE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:   W3FS15       UPDATING OFFICE NOTE 85 DATE/TIME WORD
C   PRGMMR: REJONES          ORG: NMC421     DATE: 89-08-23
C
C ABSTRACT: UPDATES OR BACKDATES A FULLWORD DATE/TIME WORD (O.N. 84)
C   BY A SPECIFIED NUMBER OF HOURS.
C
C PROGRAM HISTORY LOG:
C   ??-??-??  R.ALLARD
C   87-02-19  R.E.JONES  CLEAN UP CODE
C   87-02-19  R.E.JONES  CHANGE TO MICROSOFT FORTRAN 4.10
C   89-05-12  R.E.JONES  CORRECT ORDER OF BYTES IN DATE WORD FOR PC
C   89-08-04  R.E.JONES  CLEAN UP CODE, GET RID OF ASSIGN, CORRECTION
C                        FOR MEMORY SET TO INDEFINITE.
C   89-10-25  R.E.JONES  CHANGE TO CRAY CFT77 FORTRAN
C   95-11-15  R.E.JONES  ADD SAVE STATEMENT
C
C USAGE:    CALL W3FS15 (IDATE, JTAU, NDATE)
C   INPUT ARGUMENT LIST:
C     IDATE    - PACKED BINARY DATE/TIME AS FOLLOWS:
C                BYTE 1  IS YEAR OF CENTURY  00-99
C                BYTE 2  IS MONTH            01-12
C                BYTE 3  IS DAY OF MONTH     01-31
C                BYTE 4  IS HOUR             00-23
C                SUBROUTINE TAKES ADVANTAGE OF FORTRAN ADDRESS
C                PASSING, IDATE AND NDATE MAY BE
C                A CHARACTER*1 ARRAY OF FOUR, THE LEFT 32
C                BITS OF 64 BIT INTEGER WORD. AN OFFICE NOTE 85
C                LABEL CAN BE STORED IN
C                4 INTEGER WORDS.
C                IF INTEGER THE 2ND WORD IS USED. OUTPUT
C                IS STORED IN LEFT 32 BITS. FOR A OFFICE NOTE 84
C                LABEL THE 7TH WORD IS IN THE 4TH CRAY 64 BIT
C                INTEGER, THE LEFT 32 BITS.
C     JTAU     - INTEGER  NUMBER OF HOURS TO UPDATE (IF POSITIVE)
C                OR BACKDATE (IF NEGATIVE)
C
C   OUTPUT ARGUMENT LIST:
C     NDATE    - NEW DATE/TIME WORD RETURNED IN THE
C                SAME FORMAT AS 'IDATE'. 'NDATE' AND 'IDATE' MAY
C                BE THE SAME VARIABLE.
C
C   SUBPROGRAMS CALLED:
C     LIBRARY:
C       W3LIB    - NONE
C
C   RESTRICTIONS: THIS ROUTINE IS VALID ONLY FOR THE 20TH CENTURY.
C
C   NOTES: THE FORMAT OF THE DATE/TIME WORD IS THE SAME AS THE
C     SEVENTH WORD OF THE PACKED DATA FIELD LABEL (SEE O.N. 84) AND
C     THE THIRD WORD OF A BINARY DATA SET LABEL (SEE O.N. 85).
C
C   EXIT STATES:
C     AN ERROR FOUND BY OUT OF RANGE TESTS ON THE GIVEN DATE/TIME
C     INFORMATION WILL BE INDICATED BY RETURNING A BINARY ZERO WORD
C     IN 'NDATE'.
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN
C   MACHINE:  CRAY Y-MP8/832
C
C$$$
C
      INTEGER     ITABYR(13)
      INTEGER     LPTB(13)
      INTEGER     NOLPTB(13)
C
      CHARACTER*1 IDATE(4)
      CHARACTER*1 NDATE(4)
C
      SAVE
C
      DATA  LPTB  /0000,0744,1440,2184,2904,3648,4368,5112,
     &             5856,6576,7320,8040,8784/
      DATA  NOLPTB/0000,0744,1416,2160,2880,3624,4344,5088,
     &             5832,6552,7296,8016,8760/
      DATA  ICENTY/1900/
C
C     ...WHERE ICENTY IS FOR THE 20TH CENTURY ASSUMED FOR THE GIVEN
C     ...                 YEAR WITHIN THE CENTURY
C
      IYR    = mova2i(IDATE(1))
      IMONTH = mova2i(IDATE(2))
      IDAY   = mova2i(IDATE(3))
      IHOUR  = mova2i(IDATE(4))
C
      IF (IYR    .GT. 99) GO TO 1600
      IF (IMONTH .LE.  0) GO TO 1600
      IF (IMONTH .GT. 12) GO TO 1600
      IF (IDAY   .LE.  0) GO TO 1600
      IF (IDAY   .GT. 31) GO TO 1600
      IF (IHOUR  .LT.  0) GO TO 1600
      IF (IHOUR  .GT. 24) GO TO 1600
      IF (JTAU   .NE.  0) GO TO 100
C
        NDATE(1) = IDATE(1)
        NDATE(2) = IDATE(2)
        NDATE(3) = IDATE(3)
        NDATE(4) = IDATE(4)
        RETURN
C
  100 CONTINUE
        JAHR  = IYR + ICENTY
        KABUL = 1
        GO TO 900
C
C     ...WHERE 900 IS SUBROUTINE TO INITIALIZE ITABYR
C     ...AND RETURN THRU KABUL
C
  200 CONTINUE
        IHRYR  = IHOUR + 24 * (IDAY - 1) + ITABYR(IMONTH)
        IHRYR2 = IHRYR + JTAU
C
C     ...TO TEST FOR BACKDATED INTO PREVIOUS YEAR...
C
  300 CONTINUE
        IF (IHRYR2 .LT. 0) GO TO 700
C
      DO  400  M = 2,13
        IF (IHRYR2 .LT. ITABYR(M)) GO TO 600
  400 CONTINUE
C
C     ...IF IT FALLS THRU LOOP TO HERE, IT IS INTO NEXT YEAR...
C
        JAHR   = JAHR   + 1
        IHRYR2 = IHRYR2 - ITABYR(13)
        KABUL  = 2
        GO TO 900
C
  600 CONTINUE
        MONAT  = M      - 1
        IHRMO  = IHRYR2 - ITABYR(MONAT)
        NODAYS = IHRMO  / 24
        ITAG   = NODAYS + 1
        IUHR   = IHRMO  - NODAYS * 24
        GO TO 1500
C
C     ...ALL FINISHED.  RETURN TO CALLING PROGRAM.......................
C     ...COMES TO 700 IF NEG TOTAL HRS. BACK UP INTO PREVIOUS YEAR
C
  700 CONTINUE
        JAHR  = JAHR - 1
        KABUL = 3
        GO TO 900
C
C     ...WHICH IS CALL TO INITIALIZE ITABYR AND RETURN THRU KABUL
C
  800 CONTINUE
        IHRYR2 = ITABYR(13) + IHRYR2
        GO TO 300
C
C     ...SUBROUTINE INITYR...
C     ...CALLED BY GO TO 900 AFTER ASSIGNING RETURN NO. TO KABUL...
C     ...ITABYR HAS MONTHLY ACCUMULATING TOTAL HRS REL TO BEGIN OF YR.
C     ...DEPENDS ON WHETHER JAHR IS LEAP YEAR OR NOT.
C
  900 CONTINUE
        IQUOT  = JAHR / 4
        IRMNDR = JAHR - 4 * IQUOT
        IF (IRMNDR .NE. 0) GO TO 1000
C
C     ...WAS MODULO 4, SO MOST LIKELY A LEAP YEAR,
C
        IQUOT  = JAHR / 100
        IRMNDR = JAHR - 100 * IQUOT
        IF (IRMNDR .NE. 0) GO TO 1200
C
C     ...COMES THIS WAY IF A CENTURY YEAR...
C
        IQUOT  = JAHR / 400
        IRMNDR = JAHR - 400 * IQUOT
        IF (IRMNDR .EQ. 0) GO TO 1200
C
C     ...COMES TO 1000 IF NOT A LEAP YEAR...
C
 1000 CONTINUE
      DO  1100  I = 1,13
        ITABYR(I) = NOLPTB(I)
 1100 CONTINUE
      GO TO 1400
C
C     ...COMES TO 1200 IF LEAP YEAR
C
 1200 CONTINUE
      DO  1300  I = 1,13
        ITABYR(I) = LPTB(I)
 1300 CONTINUE
C
 1400 CONTINUE
        GO TO (200,300,800) KABUL
C
 1500 CONTINUE
        JAHR     = MOD(JAHR,100)
        NDATE(1) = CHAR(JAHR)
        NDATE(2) = CHAR(MONAT)
        NDATE(3) = CHAR(ITAG)
        NDATE(4) = CHAR(IUHR)
        RETURN
C
 1600 CONTINUE
        NDATE(1) = CHAR(0)
        NDATE(2) = CHAR(0)
        NDATE(3) = CHAR(0)
        NDATE(4) = CHAR(0)
C
C     ...WHICH FLAGS AN ERROR CONDITION ...
C
      RETURN
      END
      SUBROUTINE SABORT(STRING)
      CHARACTER*(*) STRING
      WRITE(*,*) 'ABORT:',STRING
      STOP 8
      END
