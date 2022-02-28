      SUBROUTINE SPRS(TDATA,NTDTA,NLATH,NLON,NSIG,
     *  MSAT,SFILE,NSIGSAT,ISAT,BLAT,ELAT,GROSS,SIGL)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    SPRS        FORM SUPEROBS FOR INTERACTIVE SATEMS.
C   PRGMMR: DERBER           ORG: W/NMC23    DATE: 92-07-21
C
C ABSTRACT: FORM SUPEROBS FOR INTERACTIVE SATEMS.
C
C PROGRAM HISTORY LOG:
C   92-07-21  DERBER 
C
C   INPUT ARGUMENT LIST:
C     TDATA    - OBS INFO AT OBS LOCATIONS
C     NTDTA    - NUMBER OF OBS
C     NLATH    - NUMBER OF LATITUDES POLE TO EQ. ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NSIG     - NUMBER OF LAYERS ON GAUSSIAN GRID
C     ISAT     - UNIT NUMBER INPUT SAT. STATISTICS
C     NSIGSAT  - NUMBER OF SIGMA LEVELS FOR INPUT SAT. STATS
C     BLAT     - BEGINNING LATITUDE REDUCTION OF WEIGHT REGION
C     ELAT     - ENDING LATITUDE REDUCTION OF WEIGHT REGION
C     GROSS    - PARAMETER FOR GROSS ERROR TESTING
C     SIGL     - SIGMA LAYER MIDPOINT VALUES
C
C   OUTPUT ARGUMENT LIST:
C     MSAT     - NUMBER OF SATELLITE PROFILES
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION TDATA(NTDTA,6)
C-CRA          DIMENSION RR(NSIG),SATERR(NSIG,NSIG,2)
C-CRA          DIMENSION SATERRIN(NSIGSAT,NSIGSAT,2)
C-CRA          DIMENSION SCRAT(3*NSIG)
C-CRA          DIMENSION ISCRAT(3*NSIG)
C-CRA          DIMENSION NUMTEMPS(NSIG),RR2(200),VLOC(200),ALSIG(NSIG)
C-CRA          DIMENSION INDL(NSIG)
C-CRA          DIMENSION EL(NSIG*NSIG)
C-CRA          DIMENSION AERR(NSIG*NSIG)
C-CRA          DIMENSION IST(2*NLATH+1,NLON+2)
C-CRA          DIMENSION OBSERROR(NSIG,2)
C-CRA          DIMENSION SIGL(NSIG),SIGLIN(NSIGSAT)
C-CRA          DIMENSION IOLD(NSIG)
C-CRA          DIMENSION SFILE(*)
C-CRA          REAL RRS(2*NLATH+1,NLON+2,NSIG)
C-CRA          REAL SATF(2*NLATH+1,NLON+2,NSIG)
C-CRA          LOGICAL ODIAG
 
C   WARNING!!!       DIMENSION TDATA(_NTDTA_,6) IS WRONG.
C    IT MUST USE THE DIMENSION OF THE CALLING PROGRAM.
          DIMENSION TDATA(120000,6)
          DIMENSION RR(28),SATERR(28,28,2)
          DIMENSION SATERRIN(28,28,2)
          DIMENSION SCRAT(3*28)
          DIMENSION ISCRAT(3*28)
          DIMENSION NUMTEMPS(28),RR2(200),VLOC(200),ALSIG(28)
          DIMENSION INDL(28)
          DIMENSION EL(28*28)
          DIMENSION AERR(28*28)
          DIMENSION IST(2*48+1,192+2)
          DIMENSION OBSERROR(28,2)
          DIMENSION SIGL(28),SIGLIN(28)
          DIMENSION IOLD(28)
          DIMENSION SFILE(*)
          REAL RRS(2*48+1,192+2,28)
          REAL SATF(2*48+1,192+2,28)
          LOGICAL ODIAG
C--------
C-------- LOCAL SPACE
C-------
C--------
      REWIND ISAT
      READ(ISAT,50)MSIGSAT
50    FORMAT(1X,I3)
      READ(ISAT,60)SIGLIN
60    FORMAT(1X,5E15.7)
      DO KK=1,2
       DO K=1,NSIGSAT
        READ(ISAT,60)(SATERRIN(I,K,KK),I=1,NSIGSAT)
       END DO
      END DO
      CLOSE (ISAT)
C-------------FIND OLD SIGMAS CLOSEST TO NEW SIGMAS
         DO K=1,NSIG
          LMIN=0.
          DISTMIN=1.E50
          DO L=1,NSIGSAT
           DIST=ABS(SIGL(K)-SIGLIN(L))
           IF(DIST.LT.DISTMIN) THEN
            DISTMIN=DIST
            LMIN=L
           END IF
          END DO
          IOLD(K)=LMIN
         END DO
         WRITE(6,*)' READ SATCOV AND GET IN CURRENT SIGMA COORDINATE'
         WRITE(6,*)' IOLD=',IOLD
         DO LOOP=1,2
          DO K=1,NSIG
           DO L=1,NSIG
            SATERR(K,L,LOOP)=SATERRIN(IOLD(K),IOLD(L),LOOP)
           END DO
          END DO
         END DO
C-------
      OBERMAX=-1.E50
      OBERMIN=1.E50
      RESMAX=-1.E50
      RATMAX=-1.E50
      NUMGROSS=0
C-------
C-CRA                SATF=0.
C       REAL SATF(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          SATF(ITMP,1,1)=0.
          ENDDO
C-CRA                IST=0
C       DIMENSION IST(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          IST(ITMP,1)=0
          ENDDO
      NLAT=2*NLATH
      DO 183 L=1,NSIG
 183  ALSIG(L)=FLOAT(L)
C-CRA                RRS=0.
C       REAL RRS(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RRS(ITMP,1,1)=0.
          ENDDO
      LLLCNT=0
      IACNT=1
      IPCNT=0
C-CRA                AERR=0.
C       DIMENSION AERR(NSIG*NSIG)
          DO ITMP=1,NSIG*NSIG
          AERR(ITMP)=0.
          ENDDO
C-CRA                INDL=0
C       DIMENSION INDL(NSIG)
          DO ITMP=1,NSIG
          INDL(ITMP)=0
          ENDDO
      ILON=1
      ILAT=2
      ISIG=3
      ITRES=4
      ITIME=5
      ILTYPE=6
C
C  96/3/25 - PATCH TO TRANSITION VTPR DATA ENCODEDA AS TYPE 170 TO 171
C
	  DO I=1,NTDTA 
      IF (TDATA(I,ILTYPE) .EQ. 170) TDATA(I,ILTYPE) = 171
	  ENDDO
C-CRA                NUMTEMPS=0
C       DIMENSION NUMTEMPS(NSIG),RR2(200),VLOC(200),ALSIG(NSIG)
          DO ITMP=1,NSIG
          NUMTEMPS(ITMP)=0
          ENDDO
      NREC=0
      IF(NTDTA .EQ. 0)GO TO 440
      ANLON=FLOAT(NLON)
      NNSAT=0
      DO L=1,NSIG
       OBSERROR(L,1)=SQRT(SATERR(L,L,1))
       OBSERROR(L,2)=SQRT(SATERR(L,L,2))
       WRITE(6,75321)L,OBSERROR(L,1),OBSERROR(L,2)
75321  FORMAT(' SAT ERRORS FOR L=',I3,' ARE ', 2F10.2)
      END DO
      DO 160 LLLLL=161,169
      DO 160 LLLL=1,2
      LLL=(LLLL-1)*10+LLLLL
      IC=0
      ITYPE=LLL-160
      IF(ITYPE .GT. 10) ITYPE=ITYPE-10
      IF(ITYPE .LT. 5)ISFLAG=1
      IF(ITYPE .GT. 5) ITYPE=ITYPE-5
      IF(ITYPE .EQ. 1 .OR. ITYPE .EQ. 2) IND=1
      IF(ITYPE .EQ. 3 )IND=2
      DO 100 I=1,NTDTA
        IF(IC .GE. NTDTA)GO TO 100
        IC=IC+1
        RTTYPE=TDATA(IC,ILTYPE)
        ISTYPE=NINT(RTTYPE)
        IF(ISTYPE .NE. LLL)GO TO 100
C-----------------------------------GROSS ERROR TEST ADDED HERE
        OBSERRLM=OBSERROR(MIN(NSIG,MAX(1,NINT(TDATA(IC,ISIG)))),IND)
        RESIDUAL=ABS(TDATA(IC,ITRES))
        RATIO=RESIDUAL/OBSERRLM
        OBERMAX=MAX(OBERMAX,OBSERRLM)
        OBERMIN=MIN(OBERMIN,OBSERRLM)
        RESMAX=MAX(RESMAX,RESIDUAL)
        RATMAX=MAX(RATMAX,RATIO)
        IF(RATIO.GT.GROSS) THEN
         NUMGROSS=NUMGROSS+1
         GO TO 100
        END IF
        STLAT=NINT(TDATA(IC,ILAT))
        STLON=NINT(TDATA(IC,ILON))
        NLEVSV=1
        RR2(1)=TDATA(IC,ITRES)
        VLOC(1)=TDATA(IC,ISIG)
        IF(STLON .GE. ANLON+1.)STLON=STLON-ANLON
        IF(STLON .LT. 1. )STLON=STLON+ANLON
        JLAT=STLAT
        JLON=STLON
        JLAT=MAX(1,MIN(JLAT,NLAT))
        IF(JLON .GT. NLON)JLON=JLON-NLON
        IF(IST(JLAT,JLON) .NE. 0 .AND. IST(JLAT,JLON) .LT. IND)
     *     GO TO 100
        JLATP=JLAT+1
        JLONP=JLON+1
        JLATP=MIN(JLATP,NLAT)
        DO 422 LLM=1,199
          ICP=IC+1
          IF(ICP .GT. NTDTA)GO TO 423
          IF(TDATA(ICP,ILAT) .NE. TDATA(IC,ILAT) .OR.
     *      TDATA(ICP,ILON) .NE. TDATA(IC,ILON)) GO TO 423
          ITYPEP=NINT(TDATA(ICP,ILTYPE))
          IF(ITYPEP .NE. ISTYPE) GO TO 423
          IF(TDATA(ICP,ISIG) .LT. TDATA(IC,ISIG))GO TO 423
          IC=ICP
          NLEVSV=NLEVSV+1
          RR2(NLEVSV)=TDATA(IC,ITRES)
          VLOC(NLEVSV)=TDATA(IC,ISIG)
 422    CONTINUE
 423    CONTINUE
        NLEVS=0
        IF(NLEVSV .EQ. 1) THEN
          RR(1)=RR2(1)
          INDL(1)=NINT(VLOC(1))
          NUMTEMPS(INDL(1))=NUMTEMPS(INDL(1))+1
          NLEVS=1
          GO TO 130
        END IF
        LESIG=VLOC(NLEVSV)
        LESIG=MIN(LESIG,NSIG)
        LESIG=MAX(1,LESIG)
        LBSIG=IFIX(VLOC(1))+1
        LBSIG=MIN(LBSIG,NSIG)
        LBSIG=MAX(1,LBSIG)
        IBVAL=1
        DIFFL=VLOC(1)-FLOAT(LBSIG-1)
        IF(DIFFL .LT. .02 .AND.
     *     LBSIG .GE. 2) THEN
          NLEVS=NLEVS+1
          NUMTEMPS(LBSIG-1)=NUMTEMPS(LBSIG-1)+1
          INDL(NLEVS)=LBSIG-1
          RR(NLEVS)=RR2(1)
        END IF
        IF(LBSIG .GT. LESIG) GO TO 140
        DO 129 LL=LBSIG,LESIG
          DO 131 L=IBVAL+1,NLEVSV
            IF(VLOC(L) .GT. FLOAT(LL)) GO TO 132
131       CONTINUE
          PRINT *,' INTERPOLATION ERROR ', VLOC,LBSIG,LESIG
132       IBVAL=L-1
          NLEVS=NLEVS+1
          NUMTEMPS(LL)=NUMTEMPS(LL)+1
          INDL(NLEVS)=LL
          RR(NLEVS)=((VLOC(IBVAL+1)-ALSIG(LL))*RR2(IBVAL)+
     *       (ALSIG(LL)-VLOC(IBVAL))*RR2(IBVAL+1))/
     *       (VLOC(IBVAL+1)-VLOC(IBVAL))
129     CONTINUE
140     DIFFL=LESIG+1-VLOC(NLEVSV)
        IF(DIFFL .GT. 0. .AND. DIFFL .LT. .02 .AND.
     *      LESIG .LE. NSIG-1) THEN
          NLEVS=NLEVS+1
          NUMTEMPS(LESIG+1)=NUMTEMPS(LESIG+1)+1
          INDL(NLEVS)=LESIG+1
          RR(NLEVS)=RR2(NLEVSV)
        END IF
        IF(NLEVS .GT. 0)THEN
        IST(JLAT,JLON)=IND
        DO 192 L=1,NLEVS
        SATF(JLAT,JLON,INDL(L))=SATF(JLAT,JLON,INDL(L))+1.
        RRS(JLAT,JLON,INDL(L))=RRS(JLAT,JLON,INDL(L))+RR(L)
192     CONTINUE
        ELSE
        PRINT *,' NUMBER OF LEVELS EQUAL TO ZERO ',JLAT,JLON
        END IF
130     CONTINUE
100   CONTINUE
160   CONTINUE
        DO 161 I=1,2*NLATH
        DO 161 J=1,NLON
        IF(IST(I,J) .NE. 0)THEN
        L=0
        DO 162 K=1,NSIG
        IF(SATF(I,J,K) .GT. 0.)THEN
        L=L+1
        INDL(L)=K
        RR(L)=RRS(I,J,K)/SATF(I,J,K)
        END IF
 162    CONTINUE
        IF(L .EQ. 0)GO TO 333
        STLAT=I
        WSCALE=1.
        IF(STLAT .GE. ELAT)WSCALE=2.
        IF(STLAT .GT. BLAT .AND. STLAT .LT. ELAT)THEN
          WSCALE=(2.*(STLAT-BLAT)+ELAT-STLAT)/(ELAT-BLAT)
        END IF
        NLEVS=L
        IND=IST(I,J)
        DO 101 L=1,NLEVS
        DO 101 M=1,NLEVS
          AERR(L+(M-1)*NLEVS)=SATERR(INDL(L),INDL(M),IND)*WSCALE
101     CONTINUE
C-CRA            CALL  MINV(AERR,NLEVS,NLEVS,SCRAT,DET,1.E-12,0,1)
            CALL IMINV(AERR,NLEVS,DET,ISCRAT(1),ISCRAT(NLEVS+1))
        IF(DET .LE. 1.E-12) PRINT *,' DET ERROR',IC,DET
        CALL CHLML(AERR,EL,NLEVS,NLEVS,NLEVS,ODIAG)
        IF(ODIAG)THEN
          PRINT *,I,J,NLEVS,(AERR(M),M=1,NLEVS*NLEVS)
          STOP
        END IF
        IPCNT=IPCNT+1
        DO N=1,NSIG+NSIG*NSIG+3
        SFILE(N+IACNT)=0.
        END DO
C       SFILE(1+IACNT)=NLEVS+.01
        SFILE(1+IACNT)=I+.01
        SFILE(2+IACNT)=J+.01
        SFILE(3+IACNT)=INDL(1)+.01
        IOFF=INDL(1)
        IACNT=IACNT+3
        NXIG=NSIG-IOFF+1
        DO 45 N=1,NLEVS
45      INDL(N)=INDL(N)-IOFF+1
        DO 44 N=1,NLEVS
44      SFILE(INDL(N)+IACNT)=RR(N)
        IACNT=IACNT+NSIG
        DO 444 NN=1,NLEVS
        DO 444 N=1,NLEVS
444     SFILE((INDL(NN)-1)*NXIG+INDL(N)+IACNT)=EL((NN-1)*NLEVS+N)
        IACNT=IACNT+NSIG*NSIG
        NNSAT=NNSAT+1
 333    CONTINUE
      END IF
161   CONTINUE
      PRINT *,' IACNT ',IACNT
      SFILE(1)=IPCNT
      NREC=NREC+1
      LSAT=0
      DO 350 K=1,NSIG
        LSAT=LSAT+NUMTEMPS(K)
        WRITE(6,340)NUMTEMPS(K),K
340     FORMAT(' THERE ARE ',I9,' SATEM OBS AT LEVEL ',I4)
350   CONTINUE
      PRINT *,' TOTAL NUMBER OF SATEM SPROBS=',LSAT
      PRINT *, ' NUMBER OF SATELLITE PROFILES =', NNSAT
440   MSAT=NREC
      PRINT *,' MSAT= ',MSAT
       WRITE(6,*)' GROSS ERROR CHECK FOR SATTEMPS:'
       WRITE(6,*)'   OBS ERROR MAX,MIN=',OBERMAX,OBERMIN
       WRITE(6,*)'   FOR CHECK, MAX RATIO RESIDUAL/OB ERROR =',GROSS
       WRITE(6,*)'   MAX RESIDUAL=',RESMAX
       WRITE(6,*)'   MAX RATIO=',RATMAX
       WRITE(6,*)'   NUMBER OBS THAT FAILED GROSS TEST = ',NUMGROSS
      RETURN
      END
