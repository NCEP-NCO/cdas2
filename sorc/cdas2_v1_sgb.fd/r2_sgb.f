 
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
C-----------------------------------------------------------------------
      PROGRAM SGB
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM:  SGB         TRANSFORM SIGMA TO SIGMA GRIB
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 95-01-19
C
C ABSTRACT: PROGRAM TRANSFORMS SIGMA INPUT TO SIGMA GRIB OUTPUT.
C   BY DEFAULT, TEMPERATURE, SPECIFIC HUMIDITY, DIVERENCE AND VORTICITY,
C   ZONAL AND MERIDIONAL WINDS, SURFACE PRESSURE AND ITS GRADIENT,
C   AND SURFACE OROGRAPHY AND ITS GRADIENT ARE ALL OUTPUT ON THE
C   MODEL PHYSICAL GAUSSIAN GRID AND MODEL SIGMA LEVELS.
C   PARAMETERS AND GRID OUTPUT CAN BE CONTROLLED BY NAMELIST INPUT.
C
C PROGRAM HISTORY LOG:
C   95-01-19  IREDELL
C
C NAMELISTS:
C   NAMSGB:      PARAMETERS DETERMINING OUTPUT FORMAT
C     IO         NUMBER OF LONGITUDE POINTS (DEFAULT: SET BY RDSGH)
C     JO         NUMBER OF LATITUDE POINTS (DEFAULT: SET BY RDSGH)
C     NCPUS      NUMBER OF PARALLEL PROCESSES (DEFAULT: ENVIRONMENT)
C     MXBIT      MAXIMUM NUMBER OF BITS TO PACK DATA (DEFAULT: 16)
C     IDS(255)   DECIMAL SCALING OF PACKED DATA
C                (SET TO LESS THAN -128 TO SKIP PARAMETER ALTOGETHER)
C                (DEFAULT: SET BY SUBPROGRAM IDSDEF)
C     SIGMN     MINIMUM SIGMA TO OUTPUT (DEFAULT: 0.)
C     SIGMX     MAXIMUM SIGMA TO OUTPUT (DEFAULT: 1.)
C     ICEN       FORECAST CENTER IDENTIFIER (DEFAULT: 7)
C     ICEN2      FORECAST SUB-CENTER IDENTIFIER (DEFAULT: 0)
C     IGEN       MODEL GENERATING CODE (DEFAULT: FROM SIGMA FILE)
C     IPDSX      INTEGER (IPDSX(1)+1) EXTRA PDS VALUES (DEFAULT: 0)
C                (IPDSX(1) IS NUMBER OF BYTES BEYOND 40 TO FILL IN PDS;
C                 IPDSX(2:IPDSX(1)+1) ARE VALUES OF EXTRA PDS BYTES)
C
C INPUT FILES:
C   UNIT   11    SIGMA FILE(S)
C
C OUTPUT FILES:
C   UNIT   51    SIGMA GRIB1 FILE(S)
C
C SUBPROGRAMS CALLED:
C   GNCPUS       GET ENVIRONMENT NUMBER OF PARALLEL PROCESSES
C   IDSDEF       SET DEFAULTS FOR DECIMAL SCALING
C   RDSGH        READ A SIGMA FILE HEADER
C   SGB1         TRANSFORM ONE SIGMA FILE TO SIGMA GRIB
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      PARAMETER(LEVMAX=100)
C
      DIMENSION IDATE(4)
      DIMENSION SI(LEVMAX+1),SL(LEVMAX)
      DIMENSION IDS(255),IPDSX(100)
C
      NAMELIST/NAMSGB/ IDS,SIGMN,SIGMX,
     &                 ICEN,ICEN2,IGEN,IPDSX
      DATA IDS/255*0/
      DATA SIGMN/0./,SIGMX/1./
      DATA ICEN/7/,ICEN2/0/,IGEN/0/,IPDSX/100*0/
C
      CALL GNCPUS(NCPUS)
      PRINT *,'NCPUS=',NCPUS
      CALL IDSDEF(2,IDS)
C
C  IDRT=4 .. Gaussian, IDRT=0 LAT/LON
C
      IF( 192 .EQ.36.OR. 192 .EQ.72.OR. 192 .EQ.120.OR.
     1    192 .EQ.144.OR. 192 .EQ.180.OR. 192 .EQ.360) THEN
         IDRT=0
         WRITE(6,*) 'LAT/LON GRID ASSUMED'
      ELSE
         IDRT=4
         WRITE(6,*) 'GAUSSIAN GRID ASSUMED'
      ENDIF
      READ(*,NAMSGB,END=5)
5     CONTINUE
      CALL RDSGH(11,FHOUR,IDATE,SI,SL,IRET)
      CALL SGB1(IDRT,FHOUR,IDATE,SI,SL,
     &          NCPUS,IDS,SIGMN,SIGMX,ICEN,ICEN2,IGEN,IPDSX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      STOP
      END
C-----------------------------------------------------------------------
      SUBROUTINE SGB1(IDRT,FHOUR,IDATE,SI,SL,
C-CRA&                NCPUS,IDS,SIGMN,SIGMX,ICEN,ICEN2,IGEN,IPDSX)
     &                MCPUS,IDS,SIGMN,SIGMX,ICEN,ICEN2,IGEN,IPDSX)
C
      DIMENSION IDATE(4),SI(*),SL(*),IDS(255),IPDSX(*)
C
      PARAMETER(IO= 192 , JO= 94 , MXBIT=16)
      PARAMETER(JCAP= 62 ,LEVS= 28 )
      PARAMETER(NC=( 62 +1)*( 62 +2)+1,NCTOP=( 62 +1)*2)
      PARAMETER(NFLDS=6* 28 +6)	
C
      DIMENSION SLAT((JO+1)/2),CLAT((JO+1)/2),WLAT((JO+1)/2)
      DIMENSION IFAX(20),TRIG(IO*2)
      DIMENSION EPS(NC/2),EPSTOP(NCTOP/2)
      DIMENSION S(NC,NFLDS),SSTOP(NCTOP,NFLDS)
      DIMENSION FXS(2*IO+6,NFLDS),FXYS(2*IO,(JO+1)/2,NFLDS)
C
      PARAMETER(IPUT=11,IPUU=33,IPUV=34,IPUVOR=43,IPUDIV=44,IPUQ=51)
      PARAMETER(IPUP=1,IPUZ=7,IPUPX=181,IPUPY=182,IPUZX=183,IPUZY=184)
      PARAMETER(ISFC=1,ISGLEV=107)
C
      DIMENSION IPU(NFLDS),ITL(NFLDS),IL2(NFLDS)
      CHARACTER GRIB(200+IO*JO*(MXBIT+1)/8)
C
      DIMENSION MPF(255)
C
      CALL RDSS(11,SL,IDRT,CLAT,SLAT,TRIG,IFAX,EPS,EPSTOP,S,SSTOP)
C
      JFHOUR=NINT(FHOUR)
      LEN=IO*JO
      KSZ=1
      KSD=1+LEVS
      KST=1+2*LEVS
      KSQ=1+3*LEVS
      KSPSX=1+4*LEVS
      KSPSY=2+4*LEVS
      KSU=3+4*LEVS
      KSV=3+5*LEVS
      KSPS=3+6*LEVS
      KSZS=4+6*LEVS
      KSZSX=5+6*LEVS
      KSZSY=6+6*LEVS
      DO I=KSU,KSU+LEVS-1
        IPU(I)=IPUU
      ENDDO
      DO I=KSV,KSV+LEVS-1
        IPU(I)=IPUV
      ENDDO
      DO I=KSZ,KSZ+LEVS-1
        IPU(I)=IPUVOR
      ENDDO
      DO I=KSD,KSD+LEVS-1
        IPU(I)=IPUDIV
      ENDDO
      DO I=KST,KST+LEVS-1
        IPU(I)=IPUT
      ENDDO
      DO I=KSQ,KSQ+LEVS-1
        IPU(I)=IPUQ
      ENDDO
      IPU(KSPS)=IPUP
      IPU(KSZS)=IPUZ
      IPU(KSPSX)=IPUPX
      IPU(KSPSY)=IPUPY
      IPU(KSZSX)=IPUZX
      IPU(KSZSY)=IPUZY
      DO I=KSU,KSU+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      DO I=KSV,KSV+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      DO I=KSZ,KSZ+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      DO I=KSD,KSD+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      DO I=KST,KST+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      DO I=KSQ,KSQ+LEVS-1
        ITL(I)=ISGLEV
      ENDDO
      ITL(KSPS)=ISFC
      ITL(KSZS)=ISFC
      ITL(KSPSX)=ISFC
      ITL(KSPSY)=ISFC
      ITL(KSZSX)=ISFC
      ITL(KSZSY)=ISFC
      DO I=KSU,KSU+LEVS-1
        IL2(I)=NINT(SL(I-KSU+1)*1.E4)
      ENDDO
      DO I=KSV,KSV+LEVS-1
        IL2(I)=NINT(SL(I-KSV+1)*1.E4)
      ENDDO
      DO I=KSZ,KSZ+LEVS-1
        IL2(I)=NINT(SL(I-KSZ+1)*1.E4)
      ENDDO
      DO I=KSD,KSD+LEVS-1
        IL2(I)=NINT(SL(I-KSD+1)*1.E4)
      ENDDO
      DO I=KST,KST+LEVS-1
        IL2(I)=NINT(SL(I-KST+1)*1.E4)
      ENDDO
      DO I=KSQ,KSQ+LEVS-1
        IL2(I)=NINT(SL(I-KSQ+1)*1.E4)
      ENDDO
      IL2(KSPS)=0
      IL2(KSZS)=0
      IL2(KSPSX)=0
      IL2(KSPSY)=0
      IL2(KSZSX)=0
      IL2(KSZSY)=0
      IL2MN=NINT(SIGMN*1.E4)
      IL2MX=NINT(SIGMX*1.E4)
C
      CALL MPFDEF(2,MPF)
C
      DO J=1,(JO+1)/2
        CALL TRSS(TRIG,IFAX,EPS,EPSTOP,S,SSTOP,CLAT(J),SLAT(J),FXS,J)
C
        DO I=1,IO*2
          FXS(I,KSPS)=1.E3*FXS(I,KSPS)
          DO K=1,NFLDS
            FXYS(I,J,K)=FXS(I,K)
          ENDDO
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO K=1,NFLDS
        IF(IPU(K).GT.0.AND.IDS(IPU(K)).GT.-128.AND.
     &     (ITL(K).NE.ISGLEV.OR.
     &      (IL2(K).GE.IL2MN.AND.IL2(K).LE.IL2MX))) THEN
C
          CALL POLEXT(MPF(IPU(K)),IO,FXYS(1,2,K),FXYS(1+IO,2,K),
     &                FXYS(1,1,K),FXYS(1+IO,1,K))
          CALL ROWSEP(FXYS(1,1,K))
          CALL GRIBIT(FXYS(1,1,K),IBMAP,IDRT,IO,JO,MXBIT,ACOS(SLAT(1)),
     &                28,132,ICEN,IGEN,0,IPU(K),ITL(K),0,IL2(K),
     &                IDATE(4),IDATE(2),IDATE(3),IDATE(1),1,JFHOUR,0,10,
     &                0,0,ICEN2,IDS(IPU(K)),IPDSX,
     &                0.,0.,0.,0.,0.,0.,0.,0.,
     &                GRIB,LGRIB,IERR)
          IF(IERR.EQ.0) CALL WRYTE(51,LGRIB,GRIB)
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE GNCPUS(NCPUS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: GNCPUS         GETS ENVIRONMENT NUMBER OF CPUS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 94-08-19
C
C ABSTRACT: GETS AND RETURNS THE ENVIRONMENT VARIABLE NCPUS,
C   DESIGNATING THE NUMBER OF PROCESSORS OVER WHICH TO PARALLELIZE.
C
C PROGRAM HISTORY LOG:
C   94-08-19  IREDELL
C
C USAGE:    CALL GNCPUS(NCPUS)
C   OUTPUT ARGUMENTS:
C     NCPUS        INTEGER NUMBER OF CPUS
C
C SUBPROGRAMS CALLED:
C   get_environment_variable  GET ENVIRONMENT VARIABLE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      CHARACTER*8 CNCPUS
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           NCPUS=1
C-T90      IF(1.NE.1) THEN
C-CRA        call get_environment_variable('NCPUS',CNCPUS,status=iret)
C-T90      ELSE
C-T90        IRET=0
C-T90        CALL PXFGETENV('NCPUS',5,CNCPUS,LINVAL,JRET)
C-T90        IF(JRET.EQ.0) IRET=1
C-T90      ENDIF
C-CRA      IF(IRET.EQ.1) THEN
C-CRA        READ(CNCPUS,'(BN,I8)',IOSTAT=IOS) NCPUS
C-CRA        NCPUS=MAX(NCPUS,1)
C-CRA        PRINT *,'NCPUS=',NCPUS
C-CRA      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE RDSGH(NSIG,FHOUR,IDATE,SI,SL,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    RDSGH       READ SIGMA FILE HEADER RECORD
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: READS THE HEADER RECORD FROM THE SIGMA FILE.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:      CALL RDSGH(NSIG,FHOUR,IDATE,SI,SL,IRET)
C
C   INPUT ARGUMENT LIST:
C     NSIG     - INTEGER UNIT FROM WHICH TO READ HEADER
C
C   OUTPUT ARGUMENT LIST:
C     FHOUR    - REAL FORECAST HOUR
C     IDATE    - INTEGER (4) DATE
C     SI       - REAL (LEVS+1) SIGMA INTERFACES
C     SL       - REAL (LEVS) SIGMA LEVELS
C
C   INPUT FILES:
C     NSIG     - SIGMA FILE
C
C SUBPROGRAMS CALLED:
C   MAXFAC       RETURN MAXIMUM PRIME FACTOR
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      CHARACTER*32 CLABE
      DIMENSION IDATE(4)
      DIMENSION SI( 28 +1),SL( 28 )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  READ AND EXTRACT HEADER RECORD
C  READ SIGMA SPECTRAL FILE HEADER AND DETERMINE GAUSSIAN GRID
      READ(NSIG,END=91,ERR=92) CLABE
      READ(NSIG) FHOUR,IDATE,SI,SL
      PRINT *,'IDATE=',IDATE
      RETURN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  END OF FILE ENCOUNTERED
91    IRET=1
      RETURN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  I/O ERROR ENCOUNTERED
92    IRET=2
      RETURN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
      SUBROUTINE IDSDEF(IPTV,IDS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: IDSDEF         SETS DEFAULT DECIMAL SCALINGS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: SETS DECIMAL SCALINGS DEFAULTS FOR VARIOUS PARAMETERS.
C   A DECIMAL SCALING OF -3 MEANS DATA IS PACKED IN KILO-SI UNITS.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL IDSDEF(IPTV,IDS)
C   INPUT ARGUMENTS:
C     IPTV         PARAMTER TABLE VERSION (ONLY 1 OR 2 IS RECOGNIZED)
C   OUTPUT ARGUMENTS:
C     IDS          INTEGER (255) DECIMAL SCALINGS
C                  (UNKNOWN DECIMAL SCALINGS WILL NOT BE SET)
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION IDS(255)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPTV.EQ.1.OR.IPTV.EQ.2.or.IPTV.eq.132) THEN
        IDS(001)=-1     ! PRESSURE (PA)
        IDS(002)=-1     ! SEA-LEVEL PRESSURE (PA)
        IDS(003)=5      ! PRESSURE TENDENCY (PA/S)
                        !
                        !
        IDS(006)=-1     ! GEOPOTENTIAL (M2/S2)
        IDS(007)=0      ! GEOPOTENTIAL HEIGHT (M)
        IDS(008)=0      ! GEOMETRIC HEIGHT (M)
        IDS(009)=0      ! STANDARD DEVIATION OF HEIGHT (M)
                        !
        IDS(011)=1      ! TEMPERATURE (K)
        IDS(012)=1      ! VIRTUAL TEMPERATURE (K)
        IDS(013)=1      ! POTENTIAL TEMPERATURE (K)
        IDS(014)=1      ! PSEUDO-ADIABATIC POTENTIAL TEMPERATURE (K)
        IDS(015)=1      ! MAXIMUM TEMPERATURE (K)
        IDS(016)=1      ! MINIMUM TEMPERATURE (K)
        IDS(017)=1      ! DEWPOINT TEMPERATURE (K)
        IDS(018)=1      ! DEWPOINT DEPRESSION (K)
        IDS(019)=4      ! TEMPERATURE LAPSE RATE (K/M)
        IDS(020)=0      ! VISIBILITY (M)
                        ! RADAR SPECTRA 1 ()
                        ! RADAR SPECTRA 2 ()
                        ! RADAR SPECTRA 3 ()
                        !
        IDS(025)=1      ! TEMPERATURE ANOMALY (K)
        IDS(026)=-1     ! PRESSURE ANOMALY (PA)
        IDS(027)=0      ! GEOPOTENTIAL HEIGHT ANOMALY (M)
                        ! WAVE SPECTRA 1 ()
                        ! WAVE SPECTRA 2 ()
                        ! WAVE SPECTRA 3 ()
        IDS(031)=0      ! WIND DIRECTION (DEGREES)
        IDS(032)=1      ! WIND SPEED (M/S)
        IDS(033)=1      ! ZONAL WIND (M/S)
        IDS(034)=1      ! MERIDIONAL WIND (M/S)
        IDS(035)=-4     ! STREAMFUNCTION (M2/S)
        IDS(036)=-4     ! VELOCITY POTENTIAL (M2/S)
        IDS(037)=-1     ! MONTGOMERY STREAM FUNCTION (M2/S2)
        IDS(038)=8      ! SIGMA VERTICAL VELOCITY (1/S)
        IDS(039)=3      ! PRESSURE VERTICAL VELOCITY (PA/S)
        IDS(040)=4      ! GEOMETRIC VERTICAL VELOCITY (M/S)
        IDS(041)=6      ! ABSOLUTE VORTICITY (1/S)
        IDS(042)=6      ! ABSOLUTE DIVERGENCE (1/S)
        IDS(043)=6      ! RELATIVE VORTICITY (1/S)
        IDS(044)=6      ! RELATIVE DIVERGENCE (1/S)
        IDS(045)=4      ! VERTICAL U SHEAR (1/S)
        IDS(046)=4      ! VERTICAL V SHEAR (1/S)
        IDS(047)=0      ! DIRECTION OF CURRENT (DEGREES)
                        ! SPEED OF CURRENT (M/S)
                        ! U OF CURRENT (M/S)
                        ! V OF CURRENT (M/S)
        IDS(051)=4      ! SPECIFIC HUMIDITY (KG/KG)
        IDS(052)=0      ! RELATIVE HUMIDITY (PERCENT)
        IDS(053)=4      ! HUMIDITY MIXING RATIO (KG/KG)
        IDS(054)=1      ! PRECIPITABLE WATER (KG/M2)
        IDS(055)=-1     ! VAPOR PRESSURE (PA)
        IDS(056)=-1     ! SATURATION DEFICIT (PA)
        IDS(057)=1      ! EVAPORATION (KG/M2)
        IDS(058)=1      ! CLOUD ICE (KG/M2)
        IDS(059)=6      ! PRECIPITATION RATE (KG/M2/S)
        IDS(060)=0      ! THUNDERSTORM PROBABILITY (PERCENT)
        IDS(061)=1      ! TOTAL PRECIPITATION (KG/M2)
        IDS(062)=1      ! LARGE-SCALE PRECIPITATION (KG/M2)
        IDS(063)=1      ! CONVECTIVE PRECIPITATION (KG/M2)
        IDS(064)=6      ! WATER EQUIVALENT SNOWFALL RATE (KG/M2/S)
        IDS(065)=0      ! WATER EQUIVALENT OF SNOW DEPTH (KG/M2)
        IDS(066)=2      ! SNOW DEPTH (M)
                        ! MIXED-LAYER DEPTH (M)
                        ! TRANSIENT THERMOCLINE DEPTH (M)
                        ! MAIN THERMOCLINE DEPTH (M)
                        ! MAIN THERMOCLINE ANOMALY (M)
        IDS(071)=0      ! TOTAL CLOUD COVER (PERCENT)
        IDS(072)=0      ! CONVECTIVE CLOUD COVER (PERCENT)
        IDS(073)=0      ! LOW CLOUD COVER (PERCENT)
        IDS(074)=0      ! MIDDLE CLOUD COVER (PERCENT)
        IDS(075)=0      ! HIGH CLOUD COVER (PERCENT)
        IDS(076)=1      ! CLOUD WATER (KG/M2)
                        !
        IDS(078)=1      ! CONVECTIVE SNOW (KG/M2)
        IDS(079)=1      ! LARGE SCALE SNOW (KG/M2)
        IDS(080)=1      ! WATER TEMPERATURE (K)
        IDS(081)=0      ! SEA-LAND MASK ()
                        ! DEVIATION OF SEA LEVEL FROM MEAN (M)
        IDS(083)=5      ! ROUGHNESS (M)
        IDS(084)=0      ! ALBEDO (PERCENT)
        IDS(085)=1      ! SOIL TEMPERATURE (K)
        IDS(086)=0      ! SOIL WETNESS (KG/M2)
        IDS(087)=0      ! VEGETATION (PERCENT)
                        ! SALINITY (KG/KG)
        IDS(089)=4      ! DENSITY (KG/M3)
        IDS(090)=1      ! RUNOFF (KG/M2)
        IDS(091)=0      ! ICE CONCENTRATION ()
                        ! ICE THICKNESS (M)
        IDS(093)=0      ! DIRECTION OF ICE DRIFT (DEGREES)
                        ! SPEED OF ICE DRIFT (M/S)
                        ! U OF ICE DRIFT (M/S)
                        ! V OF ICE DRIFT (M/S)
                        ! ICE GROWTH (M)
                        ! ICE DIVERGENCE (1/S)
        IDS(099)=1      ! SNOW MELT (KG/M2)
                        ! SIG HEIGHT OF WAVES AND SWELL (M)
        IDS(101)=0      ! DIRECTION OF WIND WAVES (DEGREES)
                        ! SIG HEIGHT OF WIND WAVES (M)
                        ! MEAN PERIOD OF WIND WAVES (S)
        IDS(104)=0      ! DIRECTION OF SWELL WAVES (DEGREES)
                        ! SIG HEIGHT OF SWELL WAVES (M)
                        ! MEAN PERIOD OF SWELL WAVES (S)
        IDS(107)=0      ! PRIMARY WAVE DIRECTION (DEGREES)
                        ! PRIMARY WAVE MEAN PERIOD (S)
        IDS(109)=0      ! SECONDARY WAVE DIRECTION (DEGREES)
                        ! SECONDARY WAVE MEAN PERIOD (S)
        IDS(111)=0      ! NET SOLAR RADIATIVE FLUX AT SURFACE (W/M2)
        IDS(112)=0      ! NET LONGWAVE RADIATIVE FLUX AT SURFACE (W/M2)
        IDS(113)=0      ! NET SOLAR RADIATIVE FLUX AT TOP (W/M2)
        IDS(114)=0      ! NET LONGWAVE RADIATIVE FLUX AT TOP (W/M2)
        IDS(115)=0      ! NET LONGWAVE RADIATIVE FLUX (W/M2)
        IDS(116)=0      ! NET SOLAR RADIATIVE FLUX (W/M2)
        IDS(117)=0      ! TOTAL RADIATIVE FLUX (W/M2)
                        !
                        !
                        !
        IDS(121)=0      ! LATENT HEAT FLUX (W/M2)
        IDS(122)=0      ! SENSIBLE HEAT FLUX (W/M2)
        IDS(123)=0      ! BOUNDARY LAYER DISSIPATION (W/M2)
        IDS(124)=3      ! U WIND STRESS (N/M2)
        IDS(125)=3      ! V WIND STRESS (N/M2)
                        ! WIND MIXING ENERGY (J)
                        ! IMAGE DATA ()
        IDS(128)=-1     ! MEAN SEA-LEVEL PRESSURE (STDATM) (PA)
        IDS(129)=-1     ! MEAN SEA-LEVEL PRESSURE (MAPS) (PA)
        IDS(130)=-1     ! MEAN SEA-LEVEL PRESSURE (ETA) (PA)
        IDS(131)=1      ! SURFACE LIFTED INDEX (K)
        IDS(132)=1      ! BEST LIFTED INDEX (K)
        IDS(133)=1      ! K INDEX (K)
        IDS(134)=1      ! SWEAT INDEX (K)
        IDS(135)=10     ! HORIZONTAL MOISTURE DIVERGENCE (KG/KG/S)
        IDS(136)=4      ! SPEED SHEAR (1/S)
        IDS(137)=5      ! 3-HR PRESSURE TENDENCY (PA/S)
        IDS(138)=6      ! BRUNT-VAISALA FREQUENCY SQUARED (1/S2)
        IDS(139)=11     ! POTENTIAL VORTICITY (MASS-WEIGHTED) (1/S/M)
        IDS(140)=0      ! RAIN MASK ()
        IDS(141)=0      ! FREEZING RAIN MASK ()
        IDS(142)=0      ! ICE PELLETS MASK ()
        IDS(143)=0      ! SNOW MASK ()
                        !
                        !
                        !
                        !
                        !
                        !
                        ! COVARIANCE BETWEEN V AND U (M2/S2)
                        ! COVARIANCE BETWEEN U AND T (K*M/S)
                        ! COVARIANCE BETWEEN V AND T (K*M/S)
                        !
                        !
        IDS(155)=0      ! GROUND HEAT FLUX (W/M2)
        IDS(156)=0      ! CONVECTIVE INHIBITION (W/M2)
                        ! CONVECTIVE APE (J/KG)
                        ! TURBULENT KE (J/KG)
                        ! CONDENSATION PRESSURE OF LIFTED PARCEL (PA)
        IDS(160)=0      ! CLEAR SKY UPWARD SOLAR FLUX (W/M2)
        IDS(161)=0      ! CLEAR SKY DOWNWARD SOLAR FLUX (W/M2)
        IDS(162)=0      ! CLEAR SKY UPWARD LONGWAVE FLUX (W/M2)
        IDS(163)=0      ! CLEAR SKY DOWNWARD LONGWAVE FLUX (W/M2)
        IDS(164)=0      ! CLOUD FORCING NET SOLAR FLUX (W/M2)
        IDS(165)=0      ! CLOUD FORCING NET LONGWAVE FLUX (W/M2)
        IDS(166)=0      ! VISIBLE BEAM DOWNWARD SOLAR FLUX (W/M2)
        IDS(167)=0      ! VISIBLE DIFFUSE DOWNWARD SOLAR FLUX (W/M2)
        IDS(168)=0      ! NEAR IR BEAM DOWNWARD SOLAR FLUX (W/M2)
        IDS(169)=0      ! NEAR IR DIFFUSE DOWNWARD SOLAR FLUX (W/M2)
                        !
                        !
        IDS(172)=3      ! MOMENTUM FLUX (N/M2)
        IDS(173)=0      ! MASS POINT MODEL SURFACE ()
        IDS(174)=0      ! VELOCITY POINT MODEL SURFACE ()
        IDS(175)=0      ! SIGMA LAYER NUMBER ()
        IDS(176)=2      ! LATITUDE (DEGREES)
        IDS(177)=2      ! EAST LONGITUDE (DEGREES)
                        !
                        !
                        !
        IDS(181)=9      ! X-GRADIENT LOG PRESSURE (1/M)
        IDS(182)=9      ! Y-GRADIENT LOG PRESSURE (1/M)
        IDS(183)=5      ! X-GRADIENT HEIGHT (M/M)
        IDS(184)=5      ! Y-GRADIENT HEIGHT (M/M)
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
        IDS(201)=0      ! ICE-FREE WATER SURCACE (PERCENT)
                        !
                        !
        IDS(204)=0      ! DOWNWARD SOLAR RADIATIVE FLUX (W/M2)
        IDS(205)=0      ! DOWNWARD LONGWAVE RADIATIVE FLUX (W/M2)
                        !
        IDS(207)=0      ! MOISTURE AVAILABILITY (PERCENT)
                        ! EXCHANGE COEFFICIENT (KG/M2/S)
        IDS(209)=0      ! NUMBER OF MIXED LAYER NEXT TO SFC ()
                        !
        IDS(211)=0      ! UPWARD SOLAR RADIATIVE FLUX (W/M2)
        IDS(212)=0      ! UPWARD LONGWAVE RADIATIVE FLUX (W/M2)
        IDS(213)=0      ! NON-CONVECTIVE CLOUD COVER (PERCENT)
        IDS(214)=6      ! CONVECTIVE PRECIPITATION RATE (KG/M2/S)
        IDS(215)=7      ! TOTAL DIABATIC HEATING RATE (K/S)
        IDS(216)=7      ! TOTAL RADIATIVE HEATING RATE (K/S)
        IDS(217)=7      ! TOTAL DIABATIC NONRADIATIVE HEATING RATE (K/S)
        IDS(218)=2      ! PRECIPITATION INDEX (FRACTION)
        IDS(219)=1      ! STD DEV OF IR T OVER 1X1 DEG AREA (K)
        IDS(220)=4      ! NATURAL LOG OF SURFACE PRESSURE OVER 1 KPA ()
                        !
        IDS(222)=0      ! 5-WAVE GEOPOTENTIAL HEIGHT (M)
        IDS(223)=1      ! PLANT CANOPY SURFACE WATER (KG/M2)
                        !
                        !
                        ! BLACKADARS MIXING LENGTH (M)
                        ! ASYMPTOTIC MIXING LENGTH (M)
        IDS(228)=1      ! POTENTIAL EVAPORATION (KG/M2)
        IDS(229)=0      ! SNOW PHASE-CHANGE HEAT FLUX (W/M2)
                        !
        IDS(231)=3      ! CONVECTIVE CLOUD MASS FLUX (PA/S)
        IDS(232)=0      ! DOWNWARD TOTAL RADIATION FLUX (W/M2)
        IDS(233)=0      ! UPWARD TOTAL RADIATION FLUX (W/M2)
        IDS(224)=1      ! BASEFLOW-GROUNDWATER RUNOFF (KG/M2)
        IDS(225)=1      ! STORM SURFACE RUNOFF (KG/M2)
                        !
                        !
        IDS(238)=0      ! SNOW COVER (PERCENT)
        IDS(239)=1      ! SNOW TEMPERATURE (K)
                        !
        IDS(241)=7      ! LARGE SCALE CONDENSATION HEATING RATE (K/S)
        IDS(242)=7      ! DEEP CONVECTIVE HEATING RATE (K/S)
        IDS(243)=10     ! DEEP CONVECTIVE MOISTENING RATE (KG/KG/S)
        IDS(244)=7      ! SHALLOW CONVECTIVE HEATING RATE (K/S)
        IDS(245)=10     ! SHALLOW CONVECTIVE MOISTENING RATE (KG/KG/S)
        IDS(246)=7      ! VERTICAL DIFFUSION HEATING RATE (KG/KG/S)
        IDS(247)=7      ! VERTICAL DIFFUSION ZONAL ACCELERATION (M/S/S)
        IDS(248)=7      ! VERTICAL DIFFUSION MERID ACCELERATION (M/S/S)
        IDS(249)=10     ! VERTICAL DIFFUSION MOISTENING RATE (KG/KG/S)
        IDS(250)=7      ! SOLAR RADIATIVE HEATING RATE (K/S)
        IDS(251)=7      ! LONGWAVE RADIATIVE HEATING RATE (K/S)
                        ! DRAG COEFFICIENT ()
                        ! FRICTION VELOCITY (M/S)
                        ! RICHARDSON NUMBER ()
                        !
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE RDSS(NSS,SL,IDRT,
     1                CLAT,SLAT,TRIG,IFAX,EPS,EPSTOP,SS,SSTOP)
C
      PARAMETER(IO2=2* 192 , IO22=2* 192 +6,JOHF=( 94 +1)/2)
      PARAMETER(JCAP= 62 ,LEVS= 28 )
      PARAMETER(NC=( 62 +1)*( 62 +2)+1,NCTOP=( 62 +1)*2)
      PARAMETER(NFLDS=6* 28 +6)
C
      DIMENSION SL(LEVS),CLAT(JOHF),SLAT(JOHF),TRIG(IO2),IFAX(20)
      REAL EPS((JCAP+1)*(JCAP+2)/2),EPSTOP(JCAP+1)
      REAL SS(NC,NFLDS),SSTOP(NCTOP,NFLDS)
C
      REAL WLAT(JOHF)
      REAL ENN1((JCAP+1)*(JCAP+2)),ELONN1((JCAP+1)*(JCAP+2)/2)
      REAL EON((JCAP+1)*(JCAP+2)/2),EONTOP(JCAP+1)
C
      PARAMETER(G= 9.8000E+0 ,RD= 2.8705E+2 )
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    RDSS        READ DATA FROM A SIGMA SPECTRAL FILE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: READS THE RECORDS OF OROGRAPHY, SURFACE PRESSURE,
C           DIVERGENCE AND VORTICITY, TEMPERATURE AND HUMIDITY
C           FROM A SIGMA SPECTRAL FILE.  IT IS ASSUMED THAT THE FIRST
C           TWO HEADER RECORDS OF THE FILE HAVE ALREADY BEEN READ.
C           THE GRADIENTS OF OROGRAPHY AND LOG SURFACE PRESSURE
C           AND THE WIND COMPONENTS ARE ALSO COMPUTED IN SPECTRAL SPACE.
C           THE GEOPOTENTIAL OF THE PRESSURE GRADIENT IS COMPUTED TOO.
C           ALSO, SOME SPECTRAL TRANSFORM UTILITY FIELDS ARE COMPUTED.
C           SUBPROGRAM TRSS SHOULD BE USED TO TRANSFORM TO GRID
C           AS WELL AS COMPUTE DRY TEMPERATURE AND SURFACE PRESSURE
C           AND WINDS AND GRADIENTS WITHOUT A COSINE LATITUDE FACTOR.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL RDSS(NSS,JCAP,NC,NCTOP,JOHF,IO2,LEVS,SL,
C   &                 CLAT,SLAT,TRIG,IFAX,EPS,EPSTOP,SS,SSTOP)
C
C   INPUT ARGUMENT LIST:
C     NSS      - INTEGER UNIT FROM WHICH TO READ FILE
C     JCAP     - INTEGER SPECTRAL TRUNCATION
C     NC       - INTEGER NUMBER OF SPECTRAL COEFFICIENTS
C     NCTOP    - INTEGER NUMBER OF SPECTRAL COEFFICIENTS OVER TOP
C     JOHF    - INTEGER NUMBER OF LATITUDE PAIRS IN GAUSSIAN GRID
C     IO2    - INTEGER NUMBER OF VALID DATA POINTS PER LATITUDE PAIR
C     LEVS     - INTEGER NUMBER OF LEVELS
C     SL       - REAL (LEVS) SIGMA FULL LEVEL VALUES
C
C   OUTPUT ARGUMENT LIST:
C     CLAT     - REAL (JOHF) COSINES OF LATITUDE
C     SLAT     - REAL (JOHF) SINES OF LATITUDE
C     TRIG     - REAL (IO2) TRIGONOMETRIC QUANTITIES FOR THE FFT
C     IFAX     - INTEGER (20) FACTORS FOR THE FFT
C     EPS      - REAL ((JCAP+1)*(JCAP+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C     EPSTOP   - REAL (JCAP+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C     SS       - REAL (NC,6*LEVS+6) SPECTRAL COEFS
C     SSTOP    - REAL (NCTOP,6*LEVS+6) SPECTRAL COEFS OVER TOP
C                (:,1:LEVS)             VORTICITY
C                (:,LEVS+1:2*LEVS)      DIVERGENCE
C                (:,2*LEVS+1:3*LEVS)    TEMPERATURE
C                (:,3*LEVS+1:4*LEVS)    SPECIFIC HUMIDITY
C                (:,4*LEVS+1)           D(LNPS)/DX
C                (:,4*LEVS+2)           D(LNPS)/DY
C                (:,4*LEVS+3:5*LEVS+2)  ZONAL WIND
C                (:,5*LEVS+3:6*LEVS+2)  MERIDIONAL WIND
C                (:,6*LEVS+3)           SURFACE PRESSURE
C                (:,6*LEVS+4)           OROGRAPHY
C                (:,6*LEVS+5)           D(OROG)/DX
C                (:,6*LEVS+6)           D(OROG)/DY
C
C   INPUT FILES:
C     NSS      - SIGMA SPECTRAL FILE
C
C SUBPROGRAMS CALLED:
C   ELAT         COMPUTE LATITUDES
C   FFTFAX       COMPUTE UTILITY FIELDS FOR FFT
C   GSPC         COMPUTE UTILITY FIELDS FOR SPECTRAL TRANSFORM
C   GRADQ        COMPUTE GRADIENT IN SPECTRAL SPACE
C   DZ2UV        COMPUTE VECTOR COMPONENTS IN SPECTRAL SPACE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE UTILITY FIELDS
C
      IF(IDRT.EQ.0) THEN
        CALL ELAT(JOHF,SLAT,CLAT,WLAT)
      ELSEIF(IDRT.EQ.4) THEN
        CALL GLAT(JOHF,SLAT,CLAT,WLAT)
      ELSE
        WRITE(6,*) 'IDRT should be 0 or 4 but is ',IDRT
        CALL ABORT
      ENDIF
C
C-CRA CALL FFTFAX(IO2/2,IFAX,TRIG)
      CALL    FAX(IFAX,IO2/2,3)
      CALL FFTRIG(TRIG,IO2/2,3)
C
      CALL GSPC(JCAP,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
C
C  READ SIGMA SPECTRAL DATA
C
      NR=(JCAP+1)*(JCAP+2)
      READ(NSS) (SS(I,6*LEVS+4),I=1,NR)
      READ(NSS) (SS(I,6*LEVS+3),I=1,NR)
      DO K=1,LEVS
        READ(NSS) (SS(I,2*LEVS+K),I=1,NR)
      ENDDO
      DO K=1,LEVS
        READ(NSS) (SS(I,LEVS+K),I=1,NR)
        READ(NSS) (SS(I,K),I=1,NR)
      ENDDO
      DO K=1,LEVS
        READ(NSS) (SS(I,3*LEVS+K),I=1,NR)
      ENDDO
      DO K=1,NFLDS
        DO L=0,JCAP
          SSTOP(2*L+1,K)=0.
          SSTOP(2*L+2,K)=0.
        ENDDO
      ENDDO
C
C  COMPUTE GRADIENTS AND WINDS
C
      CALL GRADQ(JCAP,ENN1,ELONN1,EON,EONTOP,SS(1,6*LEVS+4),
     &           SS(1,6*LEVS+5),SS(1,6*LEVS+6),SSTOP(1,6*LEVS+6))
      CALL GRADQ(JCAP,ENN1,ELONN1,EON,EONTOP,SS(1,6*LEVS+3),
     &           SS(1,4*LEVS+1),SS(1,4*LEVS+2),SSTOP(1,4*LEVS+2))
      DO K=1,LEVS
        CALL DZ2UV(JCAP,ENN1,ELONN1,EON,EONTOP,SS(1,LEVS+K),SS(1,K),
     &             SS(1,4*LEVS+2+K),SS(1,5*LEVS+2+K),
     &             SSTOP(1,4*LEVS+2+K),SSTOP(1,5*LEVS+2+K))
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE TRSS(TRIG,IFAX,EPS,EPSTOP,SS,SSTOP,COSLAT,SINLAT,F,J)
C
      PARAMETER(IO2=2* 192 , IO22=2* 192 +6,JOHF=( 94 +1)/2)
      PARAMETER(JCAP= 62 ,LEVS= 28 )
      PARAMETER(NC=( 62 +1)*( 62 +2)+1,NCTOP=( 62 +1)*2)
      PARAMETER(NFLDS=6* 28 +6)
C
      DIMENSION TRIG(IO2),IFAX(20)
      REAL EPS((JCAP+1)*(JCAP+2)/2),EPSTOP(JCAP+1)
      REAL SS(NC,NFLDS),SSTOP(NCTOP,NFLDS)
      REAL F(IO22,NFLDS)
C
      REAL PLN((JCAP+1)*(JCAP+2)/2),PLNTOP(JCAP+1)
      REAL WFFT(IO22,2*NFLDS)
      PARAMETER(FV= 4.6150E+2 / 2.8705E+2 -1.)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    TRSS        TRANSFORM SPECTRAL TO GRID
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: TRANSFORMS SPECTRAL TO GRIDDED DATA ON A LATITUDE PAIR
C           AND COMPUTES DRY TEMPERATURE AND SURFACE PRESSURE
C           AND WINDS AND GRADIENTS WITHOUT A COSINE LATITUDE FACTOR.
C           SUBPROGRAM RDSS SHOULD BE CALLED ALREADY
C           TO READ SPECTRAL DATA AND INITIALIZE UTILITY FIELDS.
C           THIS SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSED SEGMENT.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL TRSS(TRIG,IFAX,EPS,EPSTOP,SS,SSTOP,COSLAT,SINLAT,F,J)
C
C   INPUT ARGUMENT LIST:
C     TRIG     - REAL (IO2) TRIGONOMETRIC QUANTITIES FOR THE FFT
C     IFAX     - INTEGER (20) FACTORS FOR THE FFT
C     EPS      - REAL ((JCAP+1)*(JCAP+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C     EPSTOP   - REAL (JCAP+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C     SS       - REAL (NC,6*LEVS+6) SPECTRAL COEFS
C     SSTOP    - REAL (NCTOP,6*LEVS+6) SPECTRAL COEFS OVER TOP
C     COSLAT   - REAL COSINE OF LATITUDE OF THE NORTHERN LATITUDE
C     SINLAT   - REAL SINE OF LATITUDE OF THE NORTHERN LATITUDE
C
C   OUTPUT ARGUMENT LIST:
C     F        - REAL (IO22,6*LEVS+6) GRIDDED DATA
C                (:,1:LEVS)             VORTICITY
C                (:,1*LEVS+1:2*LEVS)    DIVERGENCE
C                (:,2*LEVS+1:3*LEVS)    TEMPERATURE
C                (:,3*LEVS+1:4*LEVS)    SPECIFIC HUMIDITY
C                (:,4*LEVS+1)           D(LNPS)/DX
C                (:,4*LEVS+2)           D(LNPS)/DY
C                (:,4*LEVS+3:5*LEVS+2)  ZONAL WIND
C                (:,5*LEVS+3:6*LEVS+2)  MERIDIONAL WIND
C                (:,6*LEVS+3)           SURFACE PRESSURE
C                (:,6*LEVS+4)           OROGRAPHY
C                (:,6*LEVS+5)           D(OROG)/DX
C                (:,6*LEVS+6)           D(OROG)/DY
C
C SUBPROGRAMS CALLED:
C   PLEG         COMPUTE ASSOCIATED LEGENDRE POLYNOMIALS
C   PSYNTH       SYNTHESIZE FOURIER FROM SPECTRAL COEFFICIENTS
C   RFFTMLT      FAST FOURIER TRANSFORM
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
C  TRANSFORM SPECTRAL COEFFICIENTS TO FOURIER COEFFICIENTS
C
      CALL PLEG(JCAP,SINLAT,COSLAT,EPS,EPSTOP,PLN,PLNTOP)
C
      CALL PSYNTH(JCAP,IO22/2,NC,NCTOP,NFLDS,PLN,PLNTOP,SS,SSTOP,F)
C
C  TRANSFORM FOURIER COEFFICIENTS TO GRIDDED DATA
C
C-CRA CALL RFFTMLT(F,WFFT,TRIG,IFAX,1,IO22/2,IO2/2,2*NFLDS,1)
      CALL FFT99M (F,WFFT,TRIG,IFAX,1,IO22/2,IO2/2,2*NFLDS,1)
C
C  MOVE SOUTHERN HEMISPHERE LATITUDE AFTER NORTHERN HEMISPHERE LATITUDE
C
      DO K=1,NFLDS
        DO I=1,IO2/2
          F(IO2/2+I,K)=F(IO22/2+I,K)
        ENDDO
      ENDDO
C
C  COMPUTE DRY TEMPERATURE FROM VIRTUAL TEMPERATURE
C  AND SURFACE PRESSURE FROM LOG SURFACE PRESSURE
C  AND DIVIDE GRADIENTS AND WINDS BY COSINE OF LATITUDE.
C
      DO K=1,LEVS
        DO I=1,IO2
          F(I,2*LEVS+K)=F(I,2*LEVS+K)/(1.+FV*F(I,3*LEVS+K))
          F(I,4*LEVS+2+K)=F(I,4*LEVS+2+K)/COSLAT
          F(I,5*LEVS+2+K)=F(I,5*LEVS+2+K)/COSLAT
        ENDDO
      ENDDO
C
      DO I=1,IO2
        F(I,6*LEVS+3)=EXP(F(I,6*LEVS+3))
        F(I,4*LEVS+1)=F(I,4*LEVS+1)/COSLAT
        F(I,4*LEVS+2)=F(I,4*LEVS+2)/COSLAT
        F(I,6*LEVS+5)=F(I,6*LEVS+5)/COSLAT
        F(I,6*LEVS+6)=F(I,6*LEVS+6)/COSLAT
      ENDDO
C
      RETURN
      END
      SUBROUTINE ROWSEP(A)
      DIMENSION A(2* 192 ,( 94 +1)/2)
      DIMENSION B( 192 , 94 )
C
      DO J=1,( 94 +1)/2
        DO I=1, 192
          B(I,J)=A(I,J)
          B(I, 94 -J+1)=A( 192 +I,J)
        ENDDO
      ENDDO
      DO IJ=1, 192 * 94
        A(IJ,1)=B(IJ,1)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
CFPP$ NOCONCUR R
      SUBROUTINE GRIBIT(F,LBM,IDRT,IM,JM,MXBIT,COLAT1,
     &                  ILPDS,IPTV,ICEN,IGEN,IBMS,IPU,ITL,IL1,IL2,
     &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,
     &                  INA,INM,ICEN2,IDS,IENS,
     &                  XLAT1,XLON1,XLAT2,XLON2,DELX,DELY,ORITRU,PROJ,
     &                  GRIB,LGRIB,IERR)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GRIBIT      CREATE GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: CREATE A GRIB MESSAGE FROM A FULL FIELD.
C   AT PRESENT, ONLY GLOBAL LATLON GRIDS AND GAUSSIAN GRIDS
C   AND REGIONAL POLAR PROJECTIONS ARE ALLOWED.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C   94-05-04  JUANG (FOR GSM AND RSM USE)
C
C USAGE:    CALL GRIBIT(F,LBM,IDRT,IM,JM,MXBIT,COLAT1,
C    &                  ILPDS,IPTV,ICEN,IGEN,IBMS,IPU,ITL,IL1,IL2,
C    &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,
C    &                  INA,INM,ICEN2,IDS,IENS,
C    &                  XLAT1,XLON1,DELX,DELY,ORITRU,PROJ,
C    &                  GRIB,LGRIB,IERR)
C   INPUT ARGUMENT LIST:
C     F        - REAL (IM*JM) FIELD DATA TO PACK INTO GRIB MESSAGE
C     LBM      - LOGICAL (IM*JM) BITMAP TO USE IF IBMS=1
C     IDRT     - INTEGER DATA REPRESENTATION TYPE
C                (0 FOR LATLON OR 4 FOR GAUSSIAN OR 5 FOR POLAR)
C     IM       - INTEGER LONGITUDINAL DIMENSION
C     JM       - INTEGER LATITUDINAL DIMENSION
C     MXBIT    - INTEGER MAXIMUM NUMBER OF BITS TO USE (0 FOR NO LIMIT)
C     COLAT1   - REAL FIRST COLATITUDE OF GRID IF IDRT=4 (RADIANS)
C     ILPDS    - INTEGER LENGTH OF THE PDS (USUALLY 28)
C     IPTV     - INTEGER PARAMETER TABLE VERSION (USUALLY 2)
C     ICEN     - INTEGER FORECAST CENTER (USUALLY 7)
C     IGEN     - INTEGER MODEL GENERATING CODE
C     IBMS     - INTEGER BITMAP FLAG (0 FOR NO BITMAP)
C     IPU      - INTEGER PARAMETER AND UNIT INDICATOR
C     ITL      - INTEGER TYPE OF LEVEL INDICATOR
C     IL1      - INTEGER FIRST LEVEL VALUE (0 FOR SINGLE LEVEL)
C     IL2      - INTEGER SECOND LEVEL VALUE
C     IYR      - INTEGER YEAR
C     IMO      - INTEGER MONTH
C     IDY      - INTEGER DAY
C     IHR      - INTEGER HOUR
C     IFTU     - INTEGER FORECAST TIME UNIT (1 FOR HOUR)
C     IP1      - INTEGER FIRST TIME PERIOD
C     IP2      - INTEGER SECOND TIME PERIOD (0 FOR SINGLE PERIOD)
C     ITR      - INTEGER TIME RANGE INDICATOR (10 FOR SINGLE PERIOD)
C     INA      - INTEGER NUMBER INCLUDED IN AVERAGE
C     INM      - INTEGER NUMBER MISSING FROM AVERAGE
C     ICEN2    - INTEGER FORECAST SUBCENTER
C                (USUALLY 0 BUT 1 FOR REANAL OR 2 FOR ENSEMBLE)
C     IDS      - INTEGER DECIMAL SCALING
C     IENS     - INTEGER (5) ENSEMBLE EXTENDED PDS VALUES
C                (APPLICATION,TYPE,IDENTIFICATION,PRODUCT,SMOOTHING)
C                (USED ONLY IF ICEN2=2 AND ILPDS>=45)
C     XLAT1    - REAL FIRST POINT OF REGIONAL LATITUDE (RADIANS)
C     XLON1    - REAL FIRST POINT OF REGIONAL LONGITUDE (RADIANS)
C     XLAT2    - REAL LAST  POINT OF REGIONAL LATITUDE (RADIANS)
C     XLON2    - REAL LAST  POINT OF REGIONAL LONGITUDE (RADIANS)
C     DELX     - REAL DX ON 60N FOR REGIONAL (M)
C     DELY     - REAL DY ON 60N FOR REGIONAL (M)
C     PROJ     - REAL POLAR PROJECTION FLAG 1 FOR NORTH -1 FOR SOUTH
C                     MERCATER PROJECTION 0
C     ORITRU   - REAL ORIENTATION OF REGIONAL POLAR PROJECTION OR
C                     TRUTH FOR REGIONAL MERCATER PROJECTION
C
C   OUTPUT ARGUMENT LIST:
C     GRIB     - CHARACTER (LGRIB) GRIB MESSAGE
C     LGRIB    - INTEGER LENGTH OF GRIB MESSAGE
C                (NO MORE THAN 100+ILPDS+IM*JM*(MXBIT+1)/8)
C     IERR     - INTEGER ERROR CODE (0 FOR SUCCESS)
C
C SUBPROGRAMS CALLED:
C   GTBITS     - COMPUTE NUMBER OF BITS AND ROUND DATA APPROPRIATELY
C   W3FI72     - ENGRIB DATA INTO A GRIB1 MESSAGE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      REAL F(IM*JM)
      LOGICAL LBM(IM*JM)
      CHARACTER GRIB(*)
      PARAMETER(NIBM= 192 * 94 )
      INTEGER IBM(NIBM),IPDS(100),IGDS(100),IBDS(100)
C-CRA INTEGER IBM(IM*JM*IBMS+1-IBMS),IPDS(100),IGDS(100),IBDS(100)
      REAL FR(NIBM)
C-CRA REAL FR(IM*JM)
C-CRA CHARACTER PDS(ILPDS)
      CHARACTER PDS(500)
C
      INTEGER IENS(5),KPROB(2),KCLUST(16),KMEMBR(80)
      REAL*8 XPROB(2)
C
      INTEGER*4 IENS4(5),KPROB4(2),KCLUST4(16),KMEMBR4(80)
      INTEGER*4 IBM4(NIBM)
      INTEGER*4 IPDS4(100),IGDS4(100),IBDS4(100)
      INTEGER*4 NBIT4,NF4,NFO4,LGRIB4,IERR4
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  DETERMINE GRID PARAMETERS
      PI=ACOS(-1.)
      NF=IM*JM
      IF(IDRT.EQ.0) THEN
        IF(IM.EQ.144.AND.JM.EQ.73) THEN
          IGRID=2
        ELSEIF(IM.EQ.360.AND.JM.EQ.181) THEN
          IGRID=3
        ELSE
          IGRID=255
        ENDIF
        IRESFL=128
        ISCAN=0
        LAT1=NINT(90.E3)
        LON1=0
        LATI=NINT(180.E3/(JM-1))
        LONI=NINT(360.E3/IM)
        IGDS09=-LAT1
        IGDS10=-LONI
        IGDS11=LATI
        IGDS12=LONI
        IGDS13=0
        IGDS14=0
      ELSEIF(IDRT.EQ.4) THEN
        IF(IM.EQ.192.AND.JM.EQ.94) THEN
          IGRID=98
        ELSEIF(IM.EQ.384.AND.JM.EQ.190) THEN
          IGRID=126
        ELSE
          IGRID=255
        ENDIF
        IRESFL=128
        ISCAN=0
        LAT1=NINT(90.E3-180.E3/PI*COLAT1)
        LON1=0
        LATI=JM/2
        LONI=NINT(360.E3/IM)
        IGDS09=-LAT1
        IGDS10=-LONI
        IGDS11=LATI
        IGDS12=LONI
        IGDS13=ISCAN
        IGDS14=0
      ELSEIF(IDRT.EQ.5) THEN    ! POLAR PROJECTION
        IGRID=255
        LAT1=NINT(180.E3/ACOS(-1.) * XLAT1)
        LON1=NINT(180.E3/ACOS(-1.) * XLON1)
        IRESFL=0
        IGDS09=NINT(ORITRU*1.E3)
        IGDS10=DELX
        IGDS11=DELY
        IF( NINT(PROJ).EQ.1  ) IGDS12=0         ! NORTH POLAR PROJ
        IF( NINT(PROJ).EQ.-1 ) IGDS12=128       ! SOUTH POLAT PROJ
        ISCAN=64
        IGDS13=ISCAN
        IGDS14=0
      ELSEIF(IDRT.EQ.1) THEN    ! MERCATER PROJECTION
        IGRID=255
        LAT1=NINT(180.E3/ACOS(-1.) * XLAT1)
        LON1=NINT(180.E3/ACOS(-1.) * XLON1)
        IRESFL=0
        IGDS09=NINT(180.E3/ACOS(-1.) * XLAT2)
        IGDS10=NINT(180.E3/ACOS(-1.) * XLON2)
        IGDS11=DELX
        IGDS12=DELY
        IGDS13=NINT(ORITRU*1.E3)
        ISCAN=64
        IGDS14=ISCAN
      ELSE
        IERR=40
        RETURN
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  RESET TIME RANGE PARAMETER IN CASE OF OVERFLOW
      IF(ITR.GE.2.AND.ITR.LE.5.AND.IP2.GE.256) THEN
        JP1=IP2
        JP2=0
        JTR=10
      ELSE
        JP1=IP1
        JP2=IP2
        JTR=ITR
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILL PDS PARAMETERS
      IPDS(01)=ILPDS    ! LENGTH OF PDS
      IPDS(02)=IPTV     ! PARAMETER TABLE VERSION ID
      IPDS(03)=ICEN     ! CENTER ID
      IPDS(04)=IGEN     ! GENERATING MODEL ID
      IPDS(05)=IGRID    ! GRID ID
      IPDS(06)=1        ! GDS FLAG
      IPDS(07)=IBMS     ! BMS FLAG
      IPDS(08)=IPU      ! PARAMETER UNIT ID
      IPDS(09)=ITL      ! TYPE OF LEVEL ID
      IPDS(10)=IL1      ! LEVEL 1 OR 0
      IPDS(11)=IL2      ! LEVEL 2
c	y2k
c     IPDS(12)=IYR      ! YEAR
      IPDS(12) = mod(IYR-1,100)+1
 
      IPDS(13)=IMO      ! MONTH
      IPDS(14)=IDY      ! DAY
      IPDS(15)=IHR      ! HOUR
      IPDS(16)=0        ! MINUTE
      IPDS(17)=IFTU     ! FORECAST TIME UNIT ID
      IPDS(18)=JP1      ! TIME PERIOD 1
      IPDS(19)=JP2      ! TIME PERIOD 2 OR 0
      IPDS(20)=JTR      ! TIME RANGE INDICATOR
      IPDS(21)=INA      ! NUMBER IN AVERAGE
      IPDS(22)=INM      ! NUMBER MISSING
c	y2k
c     IPDS(23)=20       ! CENTURY
      IPDS(23)=(IYR-1)/100+1
 
      IPDS(24)=ICEN2    ! FORECAST SUBCENTER
      IPDS(25)=IDS      ! DECIMAL SCALING
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILL GDS AND BDS PARAMETERS
      IGDS(01)=0        ! NUMBER OF VERTICAL COORDS
      IGDS(02)=255      ! VERTICAL COORD FLAG
      IGDS(03)=IDRT     ! DATA REPRESENTATION TYPE
      IGDS(04)=IM       ! EAST-WEST POINTS
      IGDS(05)=JM       ! NORTH-SOUTH POINTS
      IGDS(06)=LAT1     ! LATITUDE OF ORIGIN
      IGDS(07)=LON1     ! LONGITUDE OF ORIGIN
      IGDS(08)=IRESFL   ! RESOLUTION FLAG
      IGDS(09)=IGDS09   ! LATITUDE OF END OR ORIENTATION
      IGDS(10)=IGDS10   ! LONGITUDE OF END OR DX IN METER ON 60N
      IGDS(11)=IGDS11   ! LAT INCREMENT OR GAUSSIAN LATS OR DY IN METER
      IGDS(12)=IGDS12   ! LONGITUDE INCREMENT OR PROJECTION
      IGDS(13)=IGDS13   ! SCANNING MODE OR LAT OF INTERCUT ON EARTH FOR
      IGDS(14)=IGDS14   ! NOT USED OR SCANNING MODE FOR MERCATER
      IGDS(15)=0     ! NOT USED
      IGDS(16)=0     ! NOT USED
      IGDS(17)=0     ! NOT USED
      IGDS(18)=0     ! NOT USED
      IBDS(1)=0       ! BDS FLAGS
      IBDS(2)=0       ! BDS FLAGS
      IBDS(3)=0       ! BDS FLAGS
      IBDS(4)=0       ! BDS FLAGS
      IBDS(5)=0       ! BDS FLAGS
      IBDS(6)=0       ! BDS FLAGS
      IBDS(7)=0       ! BDS FLAGS
      IBDS(8)=0       ! BDS FLAGS
      IBDS(9)=0       ! BDS FLAGS
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILL BITMAP AND COUNT VALID DATA.  RESET BITMAP FLAG IF ALL VALID.
      NBM=NF
      IF(IBMS.NE.0) THEN
        NBM=0
        DO I=1,NF
          IF(LBM(I)) THEN
            IBM(I)=1
            NBM=NBM+1
          ELSE
            IBM(I)=0
          ENDIF
        ENDDO
        IF(NBM.EQ.NF) IPDS(7)=0
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ROUND DATA AND DETERMINE NUMBER OF BITS
      IF(NBM.EQ.0) THEN
        DO I=1,NF
          FR(I)=0.
        ENDDO
        NBIT=0
      ELSE
        CALL GTBITS(IPDS(7),IDS,NF,IBM,F,FR,FMIN,FMAX,NBIT)
C       WRITE(0,'("GTBITS:",4I4,4X,2I4,4X,2G16.6)')
C    &   IPU,ITL,IL1,IL2,IDS,NBIT,FMIN,FMAX
        IF(MXBIT.GT.0) NBIT=MIN(NBIT,MXBIT)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CREATE PRODUCT DEFINITION SECTION
      DO I=1,100
        IPDS4(I)=IPDS(I)
      ENDDO
      CALL W3FI68(IPDS4,PDS)
C-CRA CALL W3FI68(IPDS,PDS)
      IF(ICEN2.EQ.2.AND.ILPDS.GE.45) THEN
        ILAST=45
        DO I=1,5
          IENS4(I)=IENS(I)
        ENDDO
        DO I=1,2
          KPROB4(I)=KPROB(I)
        ENDDO
        DO I=1,16
          KCLUST4(I)=KCLUST(I)
        ENDDO
        DO I=1,80
          KMEMBR4(I)=KMEMBR(I)
        ENDDO
        ILAST4=ILAST
        CALL PDSENS(IENS4,KPROB4,XPROB,KCLUST4,KMEMBR4,ILAST4,PDS)
C-CRA   CALL PDSENS(IENS,KPROB,XPROB,KCLUST,KMEMBR,ILAST,PDS)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CREATE GRIB MESSAGE
      NBIT4=NBIT
      DO I=1,100
        IPDS4(I)=IPDS(I)
      ENDDO
      DO I=1,100
        IGDS4(I)=IGDS(I)
      ENDDO
      DO I=1,NF
        IBM4(I)=IBM(I)
      ENDDO
      NF4=NF
      DO I=1,100
        IBDS4(I)=IBDS(I)
      ENDDO
      CALL W3FI72(0,FR,0,NBIT4,1,IPDS4,PDS,
     &            1,255,IGDS4,0,0,IBM4,NF4,IBDS4,
     &            NFO4,GRIB,LGRIB4,IERR4)
      NFO=NFO4
      LGRIB=LGRIB4
      IERR=IERR4
C-CRA CALL W3FI72(0,FR,0,NBIT,1,IPDS,PDS,
C-CRA&            1,255,IGDS,0,0,IBM,NF,IBDS,
C-CRA&            NFO,GRIB,LGRIB,IERR)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
C-----------------------------------------------------------------------
CFPP$ NOCONCUR R
      SUBROUTINE GTBITS(IBM,IDS,LEN,MG,G,GROUND,GMIN,GMAX,NBIT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GTBITS      COMPUTE NUMBER OF BITS AND ROUND FIELD.
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: THE NUMBER OF BITS REQUIRED TO PACK A GIVEN FIELD
C   AT A PARTICULAR DECIMAL SCALING IS COMPUTED USING THE FIELD RANGE.
C   THE FIELD IS ROUNDED OFF TO THE DECIMAL SCALING FOR PACKING.
C   THE MINIMUM AND MAXIMUM ROUNDED FIELD VALUES ARE ALSO RETURNED.
C   GRIB BITMAP MASKING FOR VALID DATA IS OPTIONALLY USED.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL GTBITS(IBM,IDS,LEN,MG,G,GMIN,GMAX,NBIT)
C   INPUT ARGUMENT LIST:
C     IBM      - INTEGER BITMAP FLAG (=0 FOR NO BITMAP)
C     IDS      - INTEGER DECIMAL SCALING
C                (E.G. IDS=3 TO ROUND FIELD TO NEAREST MILLI-VALUE)
C     LEN      - INTEGER LENGTH OF THE FIELD AND BITMAP
C     MG       - INTEGER (LEN) BITMAP IF IBM=1 (0 TO SKIP, 1 TO KEEP)
C     G        - REAL (LEN) FIELD
C
C   OUTPUT ARGUMENT LIST:
C     GROUND   - REAL (LEN) FIELD ROUNDED TO DECIMAL SCALING
C                (SET TO ZERO WHERE BITMAP IS 0 IF IBM=1)
C     GMIN     - REAL MINIMUM VALID ROUNDED FIELD VALUE
C     GMAX     - REAL MAXIMUM VALID ROUNDED FIELD VALUE
C     NBIT     - INTEGER NUMBER OF BITS TO PACK
C
C SUBPROGRAMS CALLED:
C   ISRCHNE  - FIND FIRST VALUE IN AN ARRAY NOT EQUAL TO TARGET VALUE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION MG(LEN),G(LEN),GROUND(LEN)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ROUND FIELD AND DETERMINE EXTREMES WHERE BITMAP IS ON
      DS=10.**IDS
      IF(IBM.EQ.0) THEN
        GROUND(1)=NINT(G(1)*DS)/DS
        GMAX=GROUND(1)
        GMIN=GROUND(1)
        DO I=2,LEN
          GROUND(I)=NINT(G(I)*DS)/DS
          GMAX=MAX(GMAX,GROUND(I))
          GMIN=MIN(GMIN,GROUND(I))
        ENDDO
      ELSE
        I1=ISRCHNE(LEN,MG,1,0)
        IF(I1.GT.0.AND.I1.LE.LEN) THEN
          DO I=1,I1-1
            GROUND(I)=0.
          ENDDO
          GROUND(I1)=NINT(G(I1)*DS)/DS
          GMAX=GROUND(I1)
          GMIN=GROUND(I1)
          DO I=I1+1,LEN
            IF(MG(I).NE.0) THEN
              GROUND(I)=NINT(G(I)*DS)/DS
              GMAX=MAX(GMAX,GROUND(I))
              GMIN=MIN(GMIN,GROUND(I))
            ELSE
              GROUND(I)=0.
            ENDIF
          ENDDO
        ELSE
          DO I=1,LEN
            GROUND(I)=0.
          ENDDO
          GMAX=0.
          GMIN=0.
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE NUMBER OF BITS
      NBIT=LOG((GMAX-GMIN)*DS+0.9)/LOG(2.)+1.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
         FUNCTION ISRCHNE(N,IX,INCX,ITARGET)
         INTEGER IX(*),ITARGET
         J=1
         ISRCHNE=0
         IF(N.LE.0) RETURN
         IF(INCX.LT.0) J=1-(N-1)*INCX
         DO I=1,N
           IF(IX(J).NE.ITARGET) THEN
             ISRCHNE=I
             RETURN
           ENDIF
           J=J+INCX
         ENDDO
         RETURN
         END
 
      SUBROUTINE DZ2UV(M,ENN1,ELONN1,EON,EONTOP,D,Z,
     &                 U,V,UTOP,VTOP)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DZ2UV       COMPUTE WINDS FROM DIVERGENCE AND VORTICITY
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES THE WIND COMPONENTS FROM DIVERGENCE AND VORTICITY
C           IN SPECTRAL SPACE. SUBPROGRAM GSPC SHOULD BE CALLED ALREADY.
C           IF L IS THE ZONAL WAVENUMBER, N IS THE TOTAL WAVENUMBER,
C           EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) AND A IS EARTH RADIUS,
C           THEN THE ZONAL WIND COMPONENT U IS COMPUTED AS
C             U(L,N)=-I*L/(N*(N+1))*A*D(L,N)
C                    +EPS(L,N+1)/(N+1)*A*Z(L,N+1)-EPS(L,N)/N*A*Z(L,N-1)
C           AND THE MERIDIONAL WIND COMPONENT V IS COMPUTED AS
C             V(L,N)=-I*L/(N*(N+1))*A*Z(L,N)
C                    -EPS(L,N+1)/(N+1)*A*D(L,N+1)+EPS(L,N)/N*A*D(L,N-1)
C           WHERE D IS DIVERGENCE AND Z IS VORTICITY.
C           EXTRA TERMS ARE COMPUTED OVER TOP OF THE SPECTRAL TRIANGLE.
C           ADVANTAGE IS TAKEN OF THE FACT THAT EPS(L,L)=0
C           IN ORDER TO VECTORIZE OVER THE ENTIRE SPECTRAL TRIANGLE.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL DZ2UV(M,ENN1,ELONN1,EON,EONTOP,D,Z,
C    &                 U,V,UTOP,VTOP)
C
C   INPUT ARGUMENT LIST:
C     M        - INTEGER SPECTRAL TRUNCATION
C     ENN1     - REAL ((M+1)*(M+2)/2) N*(N+1)/A**2
C     ELONN1   - REAL ((M+1)*(M+2)/2) L/(N*(N+1))*A
C     EON      - REAL ((M+1)*(M+2)/2) EPSILON/N*A
C     EONTOP   - REAL (M+1) EPSILON/N*A OVER TOP
C     D        - REAL ((M+1)*(M+2)) DIVERGENCE
C     Z        - REAL ((M+1)*(M+2)) VORTICITY
C
C   OUTPUT ARGUMENT LIST:
C     U        - REAL ((M+1)*(M+2)) ZONAL WIND (TIMES COSLAT)
C     V        - REAL ((M+1)*(M+2)) MERID WIND (TIMES COSLAT)
C     UTOP     - REAL (2*(M+1)) ZONAL WIND (TIMES COSLAT) OVER TOP
C     VTOP     - REAL (2*(M+1)) MERID WIND (TIMES COSLAT) OVER TOP
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      REAL ENN1((M+1)*(M+2)/2),ELONN1((M+1)*(M+2)/2)
      REAL EON((M+1)*(M+2)/2),EONTOP(M+1)
      REAL D((M+1)*(M+2)),Z((M+1)*(M+2))
      REAL U((M+1)*(M+2)),V((M+1)*(M+2)),UTOP(2*(M+1)),VTOP(2*(M+1))
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE WINDS IN THE SPECTRAL TRIANGLE
      I=1
      U(2*I-1)=EON(I+1)*Z(2*I+1)
      U(2*I)=EON(I+1)*Z(2*I+2)
      V(2*I-1)=-EON(I+1)*D(2*I+1)
      V(2*I)=-EON(I+1)*D(2*I+2)
      DO I=2,(M+1)*(M+2)/2-1
        U(2*I-1)=ELONN1(I)*D(2*I)+EON(I+1)*Z(2*I+1)-EON(I)*Z(2*I-3)
        U(2*I)=-ELONN1(I)*D(2*I-1)+EON(I+1)*Z(2*I+2)-EON(I)*Z(2*I-2)
        V(2*I-1)=ELONN1(I)*Z(2*I)-EON(I+1)*D(2*I+1)+EON(I)*D(2*I-3)
        V(2*I)=-ELONN1(I)*Z(2*I-1)-EON(I+1)*D(2*I+2)+EON(I)*D(2*I-2)
      ENDDO
      I=(M+1)*(M+2)/2
      U(2*I-1)=ELONN1(I)*D(2*I)-EON(I)*Z(2*I-3)
      U(2*I)=-ELONN1(I)*D(2*I-1)-EON(I)*Z(2*I-2)
      V(2*I-1)=ELONN1(I)*Z(2*I)+EON(I)*D(2*I-3)
      V(2*I)=-ELONN1(I)*Z(2*I-1)+EON(I)*D(2*I-2)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE WINDS OVER TOP OF THE SPECTRAL TRIANGLE
      DO L=0,M
        I=L*(2*M+1-L)/2+M+1
        UTOP(2*L+1)=-EONTOP(L+1)*Z(2*I-1)
        UTOP(2*L+2)=-EONTOP(L+1)*Z(2*I)
        VTOP(2*L+1)=EONTOP(L+1)*D(2*I-1)
        VTOP(2*L+2)=EONTOP(L+1)*D(2*I)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE ELAT(JH,SLAT,CLAT,WLAT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    ELAT        COMPUTE EQUALLY-SPACED LATITUDE FUNCTIONS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES SINES AND COSINES AND GAUSSIAN WEIGHTS
C           OF EQUALLY-SPACED LATITUDES FROM POLE TO EQUATOR.
C           THE WEIGHTS ARE COMPUTED BASED ON ELLSAESSER (JAM,1966).
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C   93-12-28  IREDELL  MODIFIED WEIGHTS BASED ON ELLSAESSER
C   96-03-01  KANAMITSU MODIFIED FOR WORKSTATION VERSION
C
C USAGE:    CALL ELAT(JH,SLAT,CLAT,WLAT)
C
C   INPUT ARGUMENT LIST:
C     JH       - INTEGER NUMBER OF LATITUDES IN A HEMISPHERE
C
C   OUTPUT ARGUMENT LIST:
C     SLAT     - REAL (JH) SINES OF LATITUDE
C     CLAT     - REAL (JH) COSINES OF LATITUDE
C     WLAT     - REAL (JH) GAUSSIAN WEIGHTS
C
C SUBPROGRAMS CALLED:
C   MINV         SOLVES FULL MATRIX PROBLEM
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION SLAT(JH),CLAT(JH),WLAT(JH)
      PARAMETER(JJH=( 94 +1)/2)
      DIMENSION AWORK(JJH,JJH+1)
C-CRA DIMENSION BWORK(JJH*2)
      DIMENSION IWORK(JJH*2)
      PARAMETER(PI=3.14159265358979)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DLAT=0.5*PI/(JH-1)
      SLAT(1)=1.
      CLAT(1)=0.
      DO J=2,JH-1
        SLAT(J)=COS((J-1)*DLAT)
        CLAT(J)=SIN((J-1)*DLAT)
      ENDDO
      SLAT(JH)=0.
      CLAT(JH)=1.
      DO JS=1,JH
        DO J=1,JH
          AWORK(JS,J)=COS(2*(JS-1)*(J-1)*DLAT)
        ENDDO
      ENDDO
C-CRA DO JS=1,JH
C-CRA   AWORK(JS,JH+1)=-1./(4*(JS-1)**2-1)
C-CRA ENDDO
C-CRA CALL MINV (AWORK,JH,JH,BWORK,DA,1.E-12,1,0)
CCAFA CALL IMINV(AWORK,JH,1.E-12,IWORK(1),IWORK(JH+1))
      CALL IMINV(AWORK,JH,DET,IWORK(1),IWORK(JH+1))
C-CRA DO J=1,JH
C-CRA   WLAT(J)=AWORK(J,JH+1)
C-CRA ENDDO
      DO J=1,JH
        WLAT(J)=0.
      ENDDO
      DO J=1,JH
        DO JJ=1,JH
          WLAT(J)=WLAT(J)+AWORK(JJ,J)*(-1./(4.*(JJ-1)**2-1))
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE GLAT(JH,SLAT,CLAT,WLAT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GLAT        COMPUTE GAUSSIAN LATITUDE FUNCTIONS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES SINES OF GAUSSIAN LATITUDE BY ITERATION.
C           THE COSINES OF GAUSSIAN LATITUDE AND GAUSSIAN WEIGHTS
C           ARE ALSO COMPUTED.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL GLAT(JH,SLAT,CLAT,WLAT)
C
C   INPUT ARGUMENT LIST:
C     JH       - INTEGER NUMBER OF GAUSSIAN LATITUDES IN A HEMISPHERE
C
C   OUTPUT ARGUMENT LIST:
C     SLAT     - REAL (JH) SINES OF (POSITIVE) GAUSSIAN LATITUDE
C     CLAT     - REAL (JH) COSINES OF GAUSSIAN LATITUDE
C     WLAT     - REAL (JH) GAUSSIAN WEIGHTS FOR THE NH
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION SLAT(JH),CLAT(JH),WLAT(JH)
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25,EPS=1.E-14)
      PARAMETER(JBZ=50)
      DIMENSION BZ(JBZ)
      PARAMETER(JJH=( 94 +1)/2)
      DIMENSION PK(JJH),PKM1(JJH)
      DATA BZ        / 2.4048255577,  5.5200781103,
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ESTIMATE LATITUDES USING BESSEL FUNCTION
      R=1./SQRT((2*JH+0.5)**2+C)
      DO J=1,MIN(JH,JBZ)
        SLAT(J)=COS(BZ(J)*R)
      ENDDO
      DO J=JBZ+1,JH
        SLAT(J)=COS((BZ(JBZ)+(J-JBZ)*PI)*R)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CONVERGE UNTIL ALL SINES OF GAUSSIAN LATITUDE ARE WITHIN EPS
      SPMAX=1.
      DO WHILE(SPMAX.GT.EPS)
        SPMAX=0.
        DO J=1,JH
          PKM1(J)=1.
          PK(J)=SLAT(J)
        ENDDO
        DO N=2,2*JH
          DO J=1,JH
            PKM2=PKM1(J)
            PKM1(J)=PK(J)
            PK(J)=((2*N-1)*SLAT(J)*PKM1(J)-(N-1)*PKM2)/N
          ENDDO
        ENDDO
        DO J=1,JH
          SP=PK(J)*(1.-SLAT(J)**2)/(2*JH*(PKM1(J)-SLAT(J)*PK(J)))
          SLAT(J)=SLAT(J)-SP
          SPMAX=MAX(SPMAX,ABS(SP))
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE COSINES AND GAUSSIAN WEIGHTS
      DO J=1,JH
        CLAT(J)=SQRT(1.-SLAT(J)**2)
        WLAT(J)=2.*(1.-SLAT(J)**2)/(2*JH*PKM1(J))**2
      ENDDO
      RETURN
      END
      SUBROUTINE GRADQ(M,ENN1,ELONN1,EON,EONTOP,Q,
     &                 QDX,QDY,QDYTOP)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GRADQ       COMPUTE GRADIENT IN SPECTRAL SPACE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES THE HORIZONTAL VECTOR GRADIENT OF A SCALAR FIELD
C           IN SPECTRAL SPACE. SUBPROGRAM GSPC SHOULD BE CALLED ALREADY.
C           IF L IS THE ZONAL WAVENUMBER, N IS THE TOTAL WAVENUMBER,
C           EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) AND A IS EARTH RADIUS,
C           THEN THE ZONAL GRADIENT OF Q(L,N) IS SIMPLY I*L/A*Q(L,N)
C           WHILE THE MERIDIONAL GRADIENT OF Q(L,N) IS COMPUTED AS
C           EPS(L,N+1)*(N+2)/A*Q(L,N+1)-EPS(L,N+1)*(N-1)/A*Q(L,N-1).
C           EXTRA TERMS ARE COMPUTED OVER TOP OF THE SPECTRAL TRIANGLE.
C           ADVANTAGE IS TAKEN OF THE FACT THAT EPS(L,L)=0
C           IN ORDER TO VECTORIZE OVER THE ENTIRE SPECTRAL TRIANGLE.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL GRADQ(M,ENN1,ELONN1,EON,EONTOP,Q,
C    &                 QDX,QDY,QDYTOP)
C
C   INPUT ARGUMENT LIST:
C     M        - INTEGER SPECTRAL TRUNCATION
C     ENN1     - REAL ((M+1)*(M+2)/2) N*(N+1)/A**2
C     ELONN1   - REAL ((M+1)*(M+2)/2) L/(N*(N+1))*A
C     EON      - REAL ((M+1)*(M+2)/2) EPSILON/N*A
C     EONTOP   - REAL (M+1) EPSILON/N*A OVER TOP
C     Q        - REAL ((M+1)*(M+2)) SCALAR FIELD
C
C   OUTPUT ARGUMENT LIST:
C     QDX      - REAL ((M+1)*(M+2)) ZONAL GRADIENT (TIMES COSLAT)
C     QDY      - REAL ((M+1)*(M+2)) MERID GRADIENT (TIMES COSLAT)
C     QDYTOP   - REAL (2*(M+1)) MERID GRADIENT (TIMES COSLAT) OVER TOP
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      REAL ENN1((M+1)*(M+2)/2),ELONN1((M+1)*(M+2)/2)
      REAL EON((M+1)*(M+2)/2),EONTOP(M+1)
      REAL Q((M+1)*(M+2))
      REAL QDX((M+1)*(M+2)),QDY((M+1)*(M+2)),QDYTOP(2*(M+1))
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TAKE ZONAL AND MERIDIONAL GRADIENTS
      I=1
      QDX(2*I-1)=0.
      QDX(2*I)=0.
      QDY(2*I-1)=EON(I+1)*ENN1(I+1)*Q(2*I+1)
      QDY(2*I)=EON(I+1)*ENN1(I+1)*Q(2*I+2)
      DO I=2,(M+1)*(M+2)/2-1
        QDX(2*I-1)=-ELONN1(I)*ENN1(I)*Q(2*I)
        QDX(2*I)=ELONN1(I)*ENN1(I)*Q(2*I-1)
        QDY(2*I-1)=EON(I+1)*ENN1(I+1)*Q(2*I+1)-EON(I)*ENN1(I-1)*Q(2*I-3)
        QDY(2*I)=EON(I+1)*ENN1(I+1)*Q(2*I+2)-EON(I)*ENN1(I-1)*Q(2*I-2)
      ENDDO
      I=(M+1)*(M+2)/2
      QDX(2*I-1)=-ELONN1(I)*ENN1(I)*Q(2*I)
      QDX(2*I)=ELONN1(I)*ENN1(I)*Q(2*I-1)
      QDY(2*I-1)=-EON(I)*ENN1(I-1)*Q(2*I-3)
      QDY(2*I)=-EON(I)*ENN1(I-1)*Q(2*I-2)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TAKE MERIDIONAL GRADIENT OVER TOP
      DO L=0,M
        I=L*(2*M+1-L)/2+M+1
        QDYTOP(2*L+1)=-EONTOP(L+1)*ENN1(I)*Q(2*I-1)
        QDYTOP(2*L+2)=-EONTOP(L+1)*ENN1(I)*Q(2*I)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE GSPC(M,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GSPC        COMPUTE UTILITY SPECTRAL FIELDS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES CONSTANT FIELDS INDEXED IN THE SPECTRAL TRIANGLE
C           IN "IBM ORDER" (ZONAL WAVENUMBER IS THE SLOWER INDEX).
C           IF L IS THE ZONAL WAVENUMBER AND N IS THE TOTAL WAVENUMBER
C           AND A IS THE EARTH RADIUS, THEN THE FIELDS RETURNED ARE:
C           (1) NORMALIZING FACTOR EPSILON=SQRT((N**2-L**2)/(4*N**2-1))
C           (2) LAPLACIAN FACTOR N*(N+1)/A**2
C           (3) ZONAL DERIVATIVE/LAPLACIAN FACTOR L/(N*(N+1))*A
C           (4) MERIDIONAL DERIVATIVE/LAPLACIAN FACTOR EPSILON/N*A
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL GSPC(M,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
C
C   INPUT ARGUMENT LIST:
C     M        - INTEGER SPECTRAL TRUNCATION
C
C   OUTPUT ARGUMENT LIST:
C     EPS      - REAL ((M+1)*(M+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C     EPSTOP   - REAL (M+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C     ENN1     - REAL ((M+1)*(M+2)/2) N*(N+1)/A**2
C     ELONN1   - REAL ((M+1)*(M+2)/2) L/(N*(N+1))*A
C     EON      - REAL ((M+1)*(M+2)/2) EPSILON/N*A
C     EONTOP   - REAL (M+1) EPSILON/N*A OVER TOP
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      REAL EPS((M+1)*(M+2)/2),EPSTOP(M+1)
      REAL ENN1((M+1)*(M+2)/2),ELONN1((M+1)*(M+2)/2)
      REAL EON((M+1)*(M+2)/2),EONTOP(M+1)
      PARAMETER(RERTH=6.3712E6,RA2=1./RERTH**2)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO L=0,M
        ILL=L*(2*M+3-L)/2+1
        EPS(ILL)=0.
        ENN1(ILL)=RA2*L*(L+1)
        ELONN1(ILL)=RERTH/(L+1)
        EON(ILL)=0.
      ENDDO
      DO L=0,M
        IS=L*(2*M+1-L)
        IP=IS/2+1
        DO N=L+1,M
          EPS(IP+N)=SQRT(FLOAT(N**2-L**2)/FLOAT(4*N**2-1))
          ENN1(IP+N)=RA2*N*(N+1)
          ELONN1(IP+N)=RERTH*L/(N*(N+1))
          EON(IP+N)=RERTH/N*EPS(IP+N)
        ENDDO
      ENDDO
      DO L=0,M
        EPSTOP(L+1)=SQRT(FLOAT((M+1)**2-L**2)/FLOAT(4*(M+1)**2-1))
        EONTOP(L+1)=RERTH/(M+1)*EPSTOP(L+1)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE PLEG(M,SLAT,CLAT,EPS,EPSTOP,PLN,PLNTOP)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    PLEG        COMPUTE LEGENDRE POLYNOMIALS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: EVALUATES THE ORTHONORMAL ASSOCIATED LEGENDRE POLYNOMIALS
C           IN THE SPECTRAL TRIANGLE AT A GIVEN LATITUDE.
C           SUBPROGRAM GSPC SHOULD BE CALLED ALREADY.
C           IF L IS THE ZONAL WAVENUMBER, N IS THE TOTAL WAVENUMBER,
C           AND EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) THEN
C           THE FOLLOWING BOOTSTRAPPING FORMULAS ARE USED:
C           PLN(0,0)=SQRT(0.5)
C           PLN(L,L)=PLN(L-1,L-1)*CLAT*SQRT(FLOAT(2*L+1)/FLOAT(2*L))
C           PLN(L,N)=(SLAT*PLN(L,N-1)-EPS(L,N-1)*PLN(L,N-2))/EPS(L,N)
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL PLEG(M,SLAT,CLAT,EPS,EPSTOP,PLN,PLNTOP)
C
C   INPUT ARGUMENT LIST:
C     M        - INTEGER SPECTRAL TRUNCATION
C     SLAT     - REAL SINE OF LATITUDE
C     CLAT     - REAL COSINE OF LATITUDE
C     EPS      - REAL ((M+1)*(M+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C     EPSTOP   - REAL (M+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C
C   OUTPUT ARGUMENT LIST:
C     PLN      - REAL ((M+1)*(M+2)/2) LEGENDRE POLYNOMIAL
C     PLNTOP   - REAL (M+1) LEGENDRE POLYNOMIAL OVER TOP
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
CFPP$ NOCONCUR R
      REAL EPS((M+1)*(M+2)/2),EPSTOP(M+1)
      REAL PLN((M+1)*(M+2)/2),PLNTOP(M+1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ITERATIVELY COMPUTE PLN(L,L) (BOTTOM HYPOTENUSE OF TRIANGLE)
      NML=0
      I=1
      PLN(I)=SQRT(0.5)
      DO L=1,M-NML
        PLNI=PLN(I)
        I=L*(2*M+3-L)/2+(NML+1)
        PLN(I)=PLNI*CLAT*SQRT(FLOAT(2*L+1)/FLOAT(2*L))
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE PLN(L,L+1) (DIAGONAL NEXT TO BOTTOM HYPOTENUSE OF TRIANGLE)
      NML=1
      DO L=0,M-NML
        I=L*(2*M+3-L)/2+(NML+1)
        PLN(I)=SLAT*PLN(I-1)/EPS(I)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE REMAINING PLN IN SPECTRAL TRIANGLE
      DO NML=2,M
        DO L=0,M-NML
          I=L*(2*M+3-L)/2+(NML+1)
          PLN(I)=(SLAT*PLN(I-1)-EPS(I-1)*PLN(I-2))/EPS(I)
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE POLYNOMIALS OVER TOP OF SPECTRAL TRIANGLE
      DO L=0,M
        NML=M+1-L
        I=L*(2*M+3-L)/2+(NML+1)
        PLNTOP(L+1)=(SLAT*PLN(I-1)-EPS(I-1)*PLN(I-2))/EPSTOP(L+1)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE WRYTE(LU,LC,C)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    WRYTE       WRITE DATA OUT BY BYTES
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: EFFICIENTLY WRITE UNFORMATTED A CHARACETER ARRAY.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL WRYTE(LU,LC,C)
C
C   INPUT ARGUMENT LIST:
C     LU       - INTEGER UNIT TO WHICH TO WRITE
C     LC       - INTEGER NUMBER OF CHARACTERS OR BYTES TO WRITE
C     C        - CHARACETER (LC) DATA TO WRITE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      CHARACTER C(LC)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(LU) C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE PSYNTH(M,IM,NC,NCTOP,KM,PLN,PLNTOP,SPC,SPCTOP,F)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    PSYNTH      SYNTHESIZE FOURIER FROM SPECTRAL
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: SYNTHESIZES FOURIER COEFFICIENTS FROM SPECTRAL COEFFICIENTS
C           FOR A LATITUDE PAIR (NORTHERN AND SOUTHERN HEMISPHERES).
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL PSYNTH(M,IM,NC,NCTOP,KM,PLN,PLNTOP,SPC,SPCTOP,F)
C
C   INPUT ARGUMENT LIST:
C     M        - INTEGER SPECTRAL TRUNCATION
C     IM       - INTEGER DIMENSION OF FOURIER COEFFICIENTS (IM>=2*(M+1))
C     NC       - INTEGER DIMENSION OF SPECTRAL COEFFICIENTS
C                (NC>=(M+1)*(M+2))
C     NCTOP    - INTEGER DIMENSION OF SPECTRAL COEFFICIENTS OVER TOP
C                (NCTOP>=2*(M+1))
C     KM       - INTEGER NUMBER OF FIELDS
C     PLN      - REAL ((M+1)*(M+2)/2) LEGENDRE POLYNOMIAL
C     PLNTOP   - REAL (M+1) LEGENDRE POLYNOMIAL OVER TOP
C     SPC      - REAL (NC,KM) SPECTRAL COEFFICIENTS
C     SPCTOP   - REAL (NCTOP,KM) SPECTRAL COEFFICIENTS OVER TOP
C
C   OUTPUT ARGUMENT LIST:
C     F        - REAL (IM,2,KM) FOURIER COEFFICIENTS FOR LATITUDE PAIR
C
C SUBPROGRAMS CALLED:
C   SGEMVX1      CRAY LIBRARY MATRIX TIMES VECTOR
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
CFPP$ NOCONCUR R
      REAL PLN((M+1)*(M+2)/2),PLNTOP(M+1)
      REAL SPC(NC,KM),SPCTOP(NCTOP,KM)
      REAL F(IM,2,KM)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INITIALIZE FOURIER COEFFICIENTS WITH TERMS OVER TOP OF THE SPECTRUM.
C  INITIALIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
      LTOPE=MOD(M+1,2)
      LTOPO=1-LTOPE
      DO K=1,KM
        DO L=LTOPE,M,2
          F(2*L+1,1,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
          F(2*L+2,1,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
          F(2*L+1,2,K)=0.
          F(2*L+2,2,K)=0.
        ENDDO
        DO L=LTOPO,M,2
          F(2*L+1,1,K)=0.
          F(2*L+2,1,K)=0.
          F(2*L+1,2,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
          F(2*L+2,2,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FOR EACH ZONAL WAVENUMBER, SYNTHESIZE TERMS OVER TOTAL WAVENUMBER.
C  SYNTHESIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
C  COMMENTED CODE REPLACED BY LIBRARY CALLS.
      DO L=0,M
        IS=L*(2*M+1-L)
        IP=IS/2+1
        DO N=L,M,2
          DO K=1,KM
            F(2*L+1,1,K)=F(2*L+1,1,K)+PLN(IP+N)*SPC(IS+2*N+1,K)
            F(2*L+2,1,K)=F(2*L+2,1,K)+PLN(IP+N)*SPC(IS+2*N+2,K)
          ENDDO
        ENDDO
C-CRA   CALL SGEMVX1(KM,(M+2-L)/2,1.,SPC(IS+2*L+1,1),NC,4,PLN(IP+L),2,
C-CRA&               1.,F(2*L+1,1,1),IM*2)
C-CRA   CALL SGEMVX1(KM,(M+2-L)/2,1.,SPC(IS+2*L+2,1),NC,4,PLN(IP+L),2,
C-CRA&               1.,F(2*L+2,1,1),IM*2)
        DO N=L+1,M,2
          DO K=1,KM
            F(2*L+1,2,K)=F(2*L+1,2,K)+PLN(IP+N)*SPC(IS+2*N+1,K)
            F(2*L+2,2,K)=F(2*L+2,2,K)+PLN(IP+N)*SPC(IS+2*N+2,K)
          ENDDO
        ENDDO
C-CRA   CALL SGEMVX1(KM,(M+1-L)/2,1.,SPC(IS+2*L+3,1),NC,4,PLN(IP+L+1),2,
C-CRA&               1.,F(2*L+1,2,1),IM*2)
C-CRA   CALL SGEMVX1(KM,(M+1-L)/2,1.,SPC(IS+2*L+4,1),NC,4,PLN(IP+L+1),2,
C-CRA&               1.,F(2*L+2,2,1),IM*2)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SEPARATE FOURIER COEFFICIENTS FROM EACH HEMISPHERE.
C  ODD POLYNOMIALS CONTRIBUTE NEGATIVELY TO THE SOUTHERN HEMISPHERE.
      DO K=1,KM
        DO L=0,M
          F1R=F(2*L+1,1,K)
          F1I=F(2*L+2,1,K)
          F(2*L+1,1,K)=F1R+F(2*L+1,2,K)
          F(2*L+2,1,K)=F1I+F(2*L+2,2,K)
          F(2*L+1,2,K)=F1R-F(2*L+1,2,K)
          F(2*L+2,2,K)=F1I-F(2*L+2,2,K)
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ZERO OUT FOURIER WAVES OUTSIDE OF SPECTRUM
      DO L2=2*M+3,IM
        DO K=1,KM
          F(L2,1,K)=0.
          F(L2,2,K)=0.
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      SUBROUTINE MPFDEF(IPTV,MPF)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: MPFDEF         SETS DEFAULT POLE VECTOR FLAGS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: SETS FIELD IDENTIFIER DEFAULTS FOR VARIOUS PARAMETERS.
C   A FLAG OF 0 MEANS SCALAR, 1 MEANS VECTOR, AND 2 MEANS FLAG.
C   THESE IDENTIFIERS ARE USED IN INTERPOLATION.
C
C PROGRAM HISTORY LOG:
C   93-10-21  IREDELL
C
C USAGE:    CALL MPFDEF(IPTV,MPF)
C   INPUT ARGUMENTS:
C     IPTV         PARAMTER TABLE VERSION (ONLY 1 OR 2 IS RECOGNIZED)
C   OUTPUT ARGUMENTS:
C     MPF          INTEGER (255) FIELD PARAMETER IDENTIFIERS
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION MPF(255)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,255
      MPF(I)=0
      ENDDO
      IF(IPTV.EQ.1.OR.IPTV.EQ.2.OR.IPTV.eq.132) THEN
        DO I=33,34
C       MPF(033:034)=1
        MPF(I)=1
        ENDDO
        DO I=49,50
C       MPF(049:050)=1
        MPF(I)=1
        ENDDO
        DO I=95,96
C       MPF(095:096)=1
        MPF(I)=1
        ENDDO
        DO I=124,125
C       MPF(124:125)=1
        MPF(I)=1
        ENDDO
        DO I=181,182
C       MPF(181:182)=1
        MPF(I)=1
        ENDDO
        DO I=183,184
C       MPF(183:184)=1
        MPF(I)=1
        ENDDO
        DO I=247,248
C       MPF(247:248)=1
        MPF(I)=1
        ENDDO
        MPF(081)=2
        MPF(091)=2
        MPF(140)=2
        MPF(141)=2
        MPF(142)=2
        MPF(143)=2
        MPF(173)=2
        MPF(174)=2
        MPF(175)=2
        MPF(209)=2
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE POLEXT(MP,IM,FNX,FSX,FN,FS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    POLEXT      EXTRAPOLATE A FIELD TO THE POLES.
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: A GLOBAL HORIZONTAL FIELD IS EXTERPOLATED TO THE POLES.
C   POLAR SCALARS ARE THE AVERAGE OF THE CLOSEST LATITUDE CIRCLE VALUES.
C   POLAR VECTOR COMPONENTS ARE TAKEN FROM THE WAVENUMBER 1 COMPONENT
C   EXTRACTED FROM THE VALUES ON THE CLOSEST LATITUDE CIRCLE.
C   POLAR FLAGS ARE COPIED FROM THE CLOSEST PRIME MERIDIAN VALUE.
C
C PROGRAM HISTORY LOG:
C   93-04-28  IREDELL
C
C USAGE:    CALL POLEXT(MP,IM,FNX,FSX,FN,FS)
C   INPUT ARGUMENT LIST:
C     MP       - INTEGER FIELD PARAMETER IDENTIFIER
C                (0 FOR SCALAR, 1 FOR VECTOR, 2 FOR FLAG)
C     IM       - INTEGER NUMBER OF LONGITUDES
C     FNX      - REAL (IM) FIELD VALUES ON THE CLOSEST LATITUDE CIRCLE
C                TO THE NORTH POLE
C     FSX      - REAL (IM) FIELD VALUES ON THE CLOSEST LATITUDE CIRCLE
C                TO THE SOUTH POLE
C
C   OUTPUT ARGUMENT LIST:
C     FN       - REAL (IM) FIELD VALUES EXTRAPOLATED TO THE NORTH POLE
C     FS       - REAL (IM) FIELD VALUES EXTRAPOLATED TO THE SOUTH POLE
C
C ATTRIBUTES:
C   LANGUAGE: ANSI FORTRAN 77
C
C$$$
      REAL FNX(IM),FSX(IM),FN(IM),FS(IM)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET POLAR VALUES FOR SCALARS OR VECTORS
      PI=ACOS(-1.)
      IF(MP.EQ.0) THEN
C  FULL SCALAR
        FNP=0.
        FSP=0.
        DO 1010 I=1,IM
          FNP=FNP+FNX(I)
          FSP=FSP+FSX(I)
1010      CONTINUE
        FNP=FNP/IM
        FSP=FSP/IM
        DO 1020 I=1,IM
          FN(I)=FNP
          FS(I)=FSP
1020    CONTINUE
      ELSEIF(MP.EQ.1) THEN
C  FULL VECTOR
        FNPC=0.
        FNPS=0.
        FSPC=0.
        FSPS=0.
        DO 1030 I=1,IM
          CI=COS(2*PI*(I-1)/IM)
          SI=SIN(2*PI*(I-1)/IM)
          FNPC=FNPC+CI*FNX(I)
          FNPS=FNPS+SI*FNX(I)
          FSPC=FSPC+CI*FSX(I)
          FSPS=FSPS+SI*FSX(I)
1030    CONTINUE
        FNPC=2*FNPC/IM
        FNPS=2*FNPS/IM
        FSPC=2*FSPC/IM
        FSPS=2*FSPS/IM
        DO 1040 I=1,IM
          CI=COS(2*PI*(I-1)/IM)
          SI=SIN(2*PI*(I-1)/IM)
          FN(I)=FNPC*CI+FNPS*SI
          FS(I)=FSPC*CI+FSPS*SI
1040    CONTINUE
      ELSEIF(MP.EQ.2) THEN
C  FULL FLAG
        DO 1050 I=1,IM
          FN(I)=FNX(1)
          FS(I)=FSX(1)
1050    CONTINUE
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE IMINV (A,N,D,L,M)
C
C     ..................................................................
C
C        ................
C
C        PURPOSE
C           INVERT A MATRIX
C
C        USAGE
C           CALL IMINV (A,N,D,L,M)
C
C        DESCRIPTION OF PARAMETERS
C           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY
C               RESULTANT INVERSE.
C           N - ORDER OF MATRIX A
C           D - RESULTANT DETERMINANT
C           L - WORK VECTOR OF LENGTH N
C           M - WORK VECTOR OF LENGTH N
C
C        REMARKS
C           MATRIX A MUST BE A GENERAL MATRIX
C
C        .............................................
C           NONE
C
C        METHOD
C           THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT
C           IS ALSO CALCULATED. A DETERMINANT OF ZERO INDICATES THAT
C           THE MATRIX IS SINGULAR.
C
C     ..................................................................
C
      DIMENSION A(N*N),L(N),M(N)
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A, D, BIGA, HOLD
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SR........ MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMEN
C        10 MUST BE CHANGED TO DABS  .
C
C        ...............................................................
C
C        SEARCH FOR LARGEST ELEMENT
C
      D=1.0 E 0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
C  10 IF (DABS(BIGA)-DABS(A(IJ))) 15,20,20
   10 IF( ABS (BIGA)- ABS (A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C        INTERCHANGE ROWS
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI) =HOLD
C
C        INTERCHANGE COLUMNS
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI) =HOLD
C
C        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
C        CONTAINED IN BIGA)
C
   45 IF(BIGA) 48,46,48
   46 D=0.0 E 0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C        REDUCE MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=A(IK)*A(KJ)+A(IJ)
   65 CONTINUE
C
C        DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C        PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C        REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.0 E 0/BIGA
   80 CONTINUE
C
C        FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  100 K=(K-1)
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI) =HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI) =HOLD
      GO TO 100
  150 RETURN
      END
