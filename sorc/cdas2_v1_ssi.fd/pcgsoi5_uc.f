      SUBROUTINE PCGSOI(INEOFS,XHAT,XHATP,PWCON,JSAT,MSAT,
     *   NITER,MITER,JITER,JCAP,NSIG,NLATH,
     *   NLON,DEL2,WGTS,PLN,QLN,RLN,TRIGS,IFAX,IN,ISATV,
     *  NTDATA,NWDATA,NPDATA,NQDATA,NPWDAT,
     *  NTRECS,NWRECS,NPRECS,NQRECS,NPWRECS,
     *  ISCRA,NBLK,ON85DT,IOANL,INEXT,INGES,RLATS,AS,RT,RU,RV,RPW,RQ,RP,
     *  NSIGSTAT,NMDSZH,JCAPSTAT,AMPDIVT,DAMPDIVT,IDIVT,
     *  NSIGDIVT,JCAPDIVT,A3,SIGL,SIGI,DSTLAST,DSTB,ISCRA3,RRM0,
     *   MLAD,ML2LM,FACTSLM,FACTVLM,
     *   LMAD,LM2ML,FACTSML,FACTVML,
     *   RUS,RVS,RTS,RVORTS,RPLONS,RPLATS,
     *   QFILE,UVFILE,TFILE,SFILE,PWFILE,PSFILE,NSPROF)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    PCGSOI      SOLVE SPECTRAL OI EQUATION
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 91-04-02
C
C ABSTRACT: SOLVE SPECTRAL OI EQUATION. AT END OF ITERATION, CONVERT
C   TO ANALYSIS UNITS
C
C PROGRAM HISTORY LOG:
C   91-04-02  PARRISH, D., DERBER,J.
C
C   INPUT ARGUMENT LIST:
C     INEOFS   - UNIT CONTAINING FORECAST ERROR COVARIANCES
C     PWCON    - CONSTANTS FOR PRECIP. WATER CALC.
C     JSAT     - UNIT CONTAINING SAT ERROR COVARIANCE MATRICES
C     MSAT     - NUMBER OF SATELLITE PROFILES
C     NITER    - MAX NUMBER OF ITERATIONS FOR CONJUGATE GRADIENT.
C     MITER    - MAX NUMBER OF OUTER ITERATIONS
C     JITER    - CURRENT OUTER ITERATION
C     JCAP     - TRIANGULAR TRUNCATION
C     NSIG     - NUMBER OF SIGMA LEVELS
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     NLON     - NUMBER OF LONGITUDES
C     AP,BP,AQR,BQR,GR - RECURSION CONSTANTS FOR SPHERICAL HARMONICS
C     DEL2     - N*(N+1)/(A**2)
C     WGTS     - GAUSSIAN INTEGRATION WEIGHTS
C     SLAT,CLAT - SIN AND COS OF GAUSSIAN LATITUDES
C     TRIGS,IFAX - USED BY FFT
C     PE0,QE0,RO0 - STARTING FUNCTIONS FOR SPHERICAL HARMONIC RECURSIONS
C     ISATV    - ARRAY OF FLAGS DETERMINING THE USE OF SAT ERR. COV.
C     NTDATA   - NUMBER OF CONVENTIONAL TEMPERATURS
C     NWDATA   - NUMBER OF CONVENTIONAL WINDS
C     NPDATA   - NUMBER OF SURFACE PRES. DATA
C     NQDATA   - NUMBER OF SPEC. HUM. DATA
C     NPWDAT   - NUMBER OF PREC. WATER DATA
C     NTRECS   - NUMBER OF RECORDS OF CONVENTIONAL TEMPERATURS
C     NWRECS   - NUMBER OF RECORDS OF CONVENTIONAL WINDS
C     NPRECS   - NUMBER OF RECORDS OF SURFACE PRES. DATA
C     NQRECS   - NUMBER OF RECORDS OF SPEC. HUM. DATA
C     NPWRECS  - NUMBER OF RECORDS OF PREC. WATER DATA
C     ISCRA    - UNIT CONTAINING CONVENTIONAL DATA INFOR.
C     NBLK     - BLOCKING FACTOR FOR UNIT ISCRA
C     ON85DT   - O.N. 85 DATE TIME GROUP
C     IOANL    - OUTPUT UNIT NUMBER
C     INGES    - 6 HR FORECAST UNIT NUMBER
C     RLATS    - GRID LAT LOCATIONS
C     AS       - BACKGROUND ERROR COEFFICIENT
C     RT,RU,RV,RPW,RQ,RP - SCRATCH GRID ARRAYS
C     NSIGSTAT - NUMBER OF SIGMA LEVELS IN STATISTICS ARRAY
C     NMDSZH   - NUMBER OF VERTICAL MODES USED IN BALANCE EQUATION
C     JCAPSTAT - SPECTRAL TRUNCATION OF STATISTICS
C     A3       - HYDROSTATIC MATRIX
C     SIGL,SIGI - SIGMA STRUCTURE
C     DSTLAST,DSTB - DIVTEND COEFS FOR LATEST STATE AND FOR BACKGROUND
C     ISCRA3   - UNIT CONTAINING GRID STUFF NEEDED FOR TAN LIN DIVT
C     RRM0     - STARTING GRADIENT OF PENALTY FOR 1ST OUTER ITERATION.
C
C   OUTPUT ARGUMENT LIST:      
C      NONE
C
C
C REMARKS:
C
C ATTRIBUTES:
C   MACHINE:  CRAY
C
C$$$
C
C-CRA          DIMENSION PWCON(NSIG)
C-CRA          DIMENSION TRIGS(NLON*2),IFAX(10)
C-CRA          DIMENSION ISATV(NSIG),AS(4)
C-CRA          DIMENSION RUS(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION RVS(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION RTS(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION RVORTS(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION RPLONS(2*NLATH+1,NLON+2)
C-CRA          DIMENSION RPLATS(2*NLATH+1,NLON+2)
C-CRA          DIMENSION MLAD(0:JCAP,0:JCAP)
C-CRA          DIMENSION ML2LM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTSLM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTVLM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION LMAD(0:JCAP,0:JCAP)
C-CRA          DIMENSION LM2ML((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTSML((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTVML((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION QFILE(17*NQDATA),UVFILE(18*NWDATA)
C-CRA          DIMENSION TFILE(17*NTDATA)
C-CRA          DIMENSION SFILE((4+(28+2)*30)*NSPROF),PWFILE(12*NPWDAT)
C-CRA          DIMENSION PSFILE(11*NPDATA)
C-CRA          DIMENSION SIGL(NSIG),SIGI(NSIG+1)
C-CRA          DIMENSION VZ(NSIG,NSIG),VD(NSIG,NSIG)
C-CRA          DIMENSION VQ(NSIG,NSIG),VH(NSIG,NSIG)
C-CRA          DIMENSION AGVZ(0:JCAP,NSIG,NMDSZH),WGVZ(0:JCAP,NMDSZH)
C-CRA          DIMENSION BVZ(0:JCAP,NSIG,NMDSZH)
C-CRA          DIMENSION IN((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION C1(2),B1(2)
C-CRA          DIMENSION A3(NSIG,NSIG)
C-CRA          DIMENSION DSMS(NSIG),DBMS(NSIG),DIFMS(NSIG)
C-CRA          REAL RLATS(2*NLATH)
C-CRA          REAL DEL2((JCAP+1)*(JCAP+2)),WGTS(NLATH*2)
C-CRA          REAL RP((2*NLATH+1)*(NLON+2))
C-CRA          REAL RPLON((2*NLATH+1)*(NLON+2))
C-CRA          REAL RPLAT((2*NLATH+1)*(NLON+2))
C-CRA          REAL RPW((2*NLATH+1)*(NLON+2))
C-CRA          REAL RQ((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL RT((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL ST((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL RU((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL RV((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL RVORT((2*NLATH+1)*(NLON+2)*NSIG)
C-CRA          REAL PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          REAL QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          REAL RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          REAL CSHAT((JCAP+1)*(JCAP+2))
C-CRA          REAL BHALF((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          REAL BHALFP((JCAP+1)*(JCAP+2))
C-CRA          REAL BDIVT(0:JCAP,NSIG)
C-CRA          REAL XHAT((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          REAL XHATP((JCAP+1)*(JCAP+2))
C-CRA          REAL PHAT((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          REAL PHATP((JCAP+1)*(JCAP+2))
C-CRA          REAL FHAT((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          REAL FHATP((JCAP+1)*(JCAP+2))
C-CRA          REAL GHAT((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          REAL GHATP((JCAP+1)*(JCAP+2))
C-CRA          REAL DENOM(2)
C-CRA          REAL PS((JCAP+1)*(JCAP+2))
C-CRA          REAL DSTLAST((JCAP+1)*(JCAP+2),NSIG)
C-CRA          REAL DSTB((JCAP+1)*(JCAP+2),NSIG)
C-CRA          REAL FACTOR((JCAP+1)*(JCAP+2))
C-CRA          REAL FACTORI((JCAP+1)*(JCAP+2))
C-CRA          REAL BALN((JCAP+1)*(JCAP+2))
C-CRA          INTEGER IDATEG(4)
C-CRA          CHARACTER*4 ON85DT(8)
C-CRA          CHARACTER*4 ON85(8)
 
          DIMENSION PWCON(28)
          DIMENSION TRIGS(192*2),IFAX(10)
          DIMENSION ISATV(28),AS(4)
          DIMENSION RUS(2*48+1,192+2,28)
          DIMENSION RVS(2*48+1,192+2,28)
          DIMENSION RTS(2*48+1,192+2,28)
          DIMENSION RVORTS(2*48+1,192+2,28)
          DIMENSION RPLONS(2*48+1,192+2)
          DIMENSION RPLATS(2*48+1,192+2)
          DIMENSION MLAD(0:62,0:62)
          DIMENSION ML2LM((62+1)*(62+2))
          DIMENSION FACTSLM((62+1)*(62+2))
          DIMENSION FACTVLM((62+1)*(62+2))
          DIMENSION LMAD(0:62,0:62)
          DIMENSION LM2ML((62+1)*(62+2))
          DIMENSION FACTSML((62+1)*(62+2))
          DIMENSION FACTVML((62+1)*(62+2))
          DIMENSION QFILE(17*15000),UVFILE(18*85000)
          DIMENSION TFILE(17*60000)
          DIMENSION SFILE((4+(28+2)*30)*10000)
          DIMENSION PWFILE(12*1),PSFILE(11*18000)
          DIMENSION SIGL(28),SIGI(28+1)
          DIMENSION VZ(28,28),VD(28,28)
          DIMENSION VQ(28,28),VH(28,28)
          DIMENSION AGVZ(0:62,28,28)
          DIMENSION WGVZ(0:62,28)
          DIMENSION BVZ(0:62,28,28)
          DIMENSION IN((62+1)*(62+2))
          DIMENSION C1(2),_B1(2)
          DIMENSION A3(28,28)
          DIMENSION DSMS(28),DBMS(28),DIFMS(28)
          REAL RLATS(2*48)
          REAL DEL2((62+1)*(62+2)),WGTS(48*2)
          REAL RP((2*48+1)*(192+2))
          REAL RPLON((2*48+1)*(192+2))
          REAL RPLAT((2*48+1)*(192+2))
          REAL RPW((2*48+1)*(192+2))
          REAL RQ((2*48+1)*(192+2)*28)
          REAL RT((2*48+1)*(192+2)*28)
          REAL ST((2*48+1)*(192+2)*28)
          REAL RU((2*48+1)*(192+2)*28)
          REAL RV((2*48+1)*(192+2)*28)
          REAL RVORT((2*48+1)*(192+2)*28)
          REAL PLN((62+1)*(62+2),48)
          REAL QLN((62+1)*(62+2),48)
          REAL RLN((62+1)*(62+2),48)
          REAL CSHAT((62+1)*(62+2))
          REAL BHALF((62+1)*(62+2),28,4)
          REAL BHALFP((62+1)*(62+2))
          REAL BDIVT(0:62,28)
          REAL XHAT((62+1)*(62+2),28,4)
          REAL XHATP((62+1)*(62+2))
          REAL PHAT((62+1)*(62+2),28,4)
          REAL PHATP((62+1)*(62+2))
          REAL FHAT((62+1)*(62+2),28,4)
          REAL FHATP((62+1)*(62+2))
          REAL GHAT((62+1)*(62+2),28,4)
          REAL GHATP((62+1)*(62+2))
          REAL DENOM(2)
          REAL PS((62+1)*(62+2))
          REAL DSTLAST((62+1)*(62+2),28)
          REAL DSTB((62+1)*(62+2),28)
          REAL FACTOR((62+1)*(62+2))
          REAL FACTORI((62+1)*(62+2))
          REAL BALN((62+1)*(62+2))
          INTEGER IDATEG(4)
          CHARACTER*4 ON85DT(8)
          CHARACTER*4 ON85(8)
C--------
C-------- LOCAL SPACE
C--------
C-------------
      NCEF=(JCAP+1)*(JCAP+2)
      NCEFS=NCEF*NSIG
      NCEF4S=4*NCEFS
      NGRP=(2*NLATH+1)*(NLON+2)
      NGRPS=(2*NLATH+1)*(NLON+2)*NSIG
C     IMAX=222201
C--------
C--------
C-------- OBTAIN INITIAL R,P,X
C--------
C--------
C-------- 7.  OBTAIN BHALF (BACKGROUND ERROR COVAR)
C--------
         II=-1
         DO M=0,JCAP
          II=II+2
          FACTOR(II)=1.
          FACTORI(II)=2.
          FACTOR(II+1)=0.
          FACTORI(II+1)=0.
          IF(M.LT.JCAP) THEN
           DO L=1,JCAP-M
            II=II+2
            FACTOR(II)=2.
            FACTORI(II)=1.
            FACTOR(II+1)=2.
            FACTORI(II+1)=1.
           END DO
          END IF
         END DO
      FACTOR(1)=0.
      FACTORI(1)=0.
C-CRA                FACTOR=FACTOR/(16.*ATAN(1.))
C       REAL FACTOR((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          FACTOR(ITMP)=FACTOR(ITMP)/(16.*ATAN(1.))
          ENDDO
      CALL GETBALN(BALN,JCAP)
C-CRA                RP=0.
C       REAL RP((2*NLATH+1)*(NLON+2))
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          RP(ITMP)=0.
          ENDDO
      CALL INITPS(RP,NLATH,NLON,NPRECS,PSFILE)
C-CRA                RU=0.
C       REAL RU((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RU(ITMP)=0.
          ENDDO
C-CRA                RV=0.
C       REAL RV((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RV(ITMP)=0.
          ENDDO
      CALL INITW(RU,RV,NLATH,NLON,NSIG,NWRECS,UVFILE)
      RLKM=400.
      CALL SATCOV(JCAP,NLATH,NLON,CSHAT,RLATS,
     *           PLN,WGTS,TRIGS,IFAX,RLKM,LMAD)
C-CRA                RT=0.
C       REAL RT((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RT(ITMP)=0.
          ENDDO
      CALL INITT(ST,NTRECS,NLATH,NLON,NSIG,TFILE)
C-CRA                RQ=0.
C       REAL RQ((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RQ(ITMP)=0.
          ENDDO
      CALL INITQPW(RQ,NLATH,NLON,NSIG,NQRECS,NPWRECS,
     *   PWCON,QFILE,PWFILE)
      CALL INITSAT(MSAT,NLATH,NLON,NSIG,RT,CSHAT,PLN,
     *  TRIGS,IFAX,JCAP,ISATV,SFILE)
C-CRA                RT=RT+ST
C       REAL RT((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RT(ITMP)=RT(ITMP)+ST(ITMP)
          ENDDO
C--------
C        PASSED THIS PLACE
C--------  FIRST OBTAIN VERTICAL RESOLUTION OF INPUT STATS
C--------
      CALL GTBHALF(INEOFS,BHALF,BHALFP,JCAP,NSIG,NLATH,AS,JCAPSTAT,
     *     NSIGSTAT,AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,SIGL)
C
C  
C
      CALL GTBDIVT(IDIVT,BDIVT,NSIG,JCAP,NSIGDIVT,JCAPDIVT,SIGL)
      PRINT *,' IN PCGSOI, AMPDIVT=',AMPDIVT
C
C     PASSED THIS PLACE
C
C-CRA                BDIVT=AMPDIVT*BDIVT
C       REAL BDIVT(0:JCAP,NSIG)
          DO JTMP=1,NSIG
          DO ITMP=0,JCAP
          BDIVT(ITMP,JTMP)=AMPDIVT*BDIVT(ITMP,JTMP)
          ENDDO
          ENDDO
      IF(JITER.GT.1) THEN
C-CRA                 DBMS=0.
C       DIMENSION DSMS(NSIG),DBMS(NSIG),DIFMS(NSIG)
          DO ITMP=1,NSIG
          DBMS(ITMP)=0.
          ENDDO
C-CRA                 DSMS=0.
C       DIMENSION DSMS(NSIG),DBMS(NSIG),DIFMS(NSIG)
          DO ITMP=1,NSIG
          DSMS(ITMP)=0.
          ENDDO
C-CRA                 DIFMS=0.
C       DIMENSION DSMS(NSIG),DBMS(NSIG),DIFMS(NSIG)
          DO ITMP=1,NSIG
          DIFMS(ITMP)=0.
          ENDDO
       DBTOT=0.
       DSTOT=0.
       DIFTOT=0.
       DO K=1,NSIG
        DO I=1,NCEF
         DSMS(K)=DSMS(K)+FACTOR(I)*DSTLAST(I,K)**2
         DBMS(K)=DBMS(K)+FACTOR(I)*DSTB(I,K)**2
         DIFMS(K)=DIFMS(K)+FACTOR(I)*(DSTLAST(I,K)-DSTB(I,K))**2
         DSTLAST(I,K)=DAMPDIVT*DSTB(I,K)-DSTLAST(I,K)
C-------------------DERBER ALWAYS CHANGES THE STUPID SIGN
         DSTLAST(I,K)=-DSTLAST(I,K)
        END DO
        DBTOT=DBTOT+DBMS(K)/NSIG
        DSTOT=DSTOT+DSMS(K)/NSIG
        DIFTOT=DIFTOT+DIFMS(K)/NSIG
       END DO
       DO K=1,NSIG
        DO I=1,NCEF
         DSTLAST(I,K)=BDIVT(IN(I),K)*DSTLAST(I,K)
        END DO
C-------------------MULTIPLY ZONAL TERMS BY EXTRA FACTOR OF 2
C------------- FOR HOMOGENEOUS, ISOTROPIC COVARIANCE IN GRID SPACE.
C-CRA                  DSTLAST(1:NCEF,K)=FACTORI(1:NCEF)*DSTLAST(1:NCEF,K)
          DO ITMP=1,NCEF
                  DSTLAST(ITMP,K)=FACTORI(ITMP)*DSTLAST(ITMP,K)
          END DO
       END DO
       WRITE(6,*)' STATS FOR DIV TEND FOLLOW FOR OUTER ITERATION =',
     *              JITER
       WRITE(6,*)' CURRENT ANALYSIS DIVTEND   BACKGROUND DIVTEND  DIFF'
       DO K=1,NSIG
        WRITE(6,94512)K,DSMS(K),DBMS(K),DIFMS(K)
94512   FORMAT(' K=',I3,3E19.5)
       END DO
       K=NSIG+1
       WRITE(6,94512)K,DSTOT,DBTOT,DIFTOT
C
C      THIS PART DID NOT PASS
C
       IF(AMPDIVT.GT.0.)
     *  CALL QTOPER(RU,RV,RVORT,RT,RPLON,RPLAT,
     *    NSIG,JCAP,NLON,NLATH,PLN,QLN,RLN,TRIGS,IFAX,
     *    DEL2,WGTS,A3,SIGL,SIGI,DSTLAST,ISCRA3,RLATS,
     *    RUS,RVS,RTS,RVORTS,RPLONS,RPLATS)
      ELSE
C-CRA                 RVORT=0.
C       REAL RVORT((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RVORT(ITMP)=0.
          ENDDO
C-CRA                 RPLAT=0.
C       REAL RPLAT((2*NLATH+1)*(NLON+2))
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          RPLAT(ITMP)=0.
          ENDDO
C-CRA                 RPLON=0.
C       REAL RPLON((2*NLATH+1)*(NLON+2))
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          RPLON(ITMP)=0.
          ENDDO
      END IF
C
      PRINT *,' BEFORE NEW HTOPER, ',RU(1),RV(1),RVORT(1),
     *   RT(1),RP(1),RPLON(1),RPLAT(1),RQ(1)
      CALL HTOPER(GHAT(1,1,1),GHAT(1,1,2),GHAT(1,1,3),
     *    GHAT(1,1,4),GHATP,RU,RV,RVORT,RT,RP,RPLON,RPLAT,RQ,
     *    BHALF,BHALFP,NSIG,JCAP,NLON,NLATH,DEL2,
     *    PLN,QLN,RLN,TRIGS,IFAX,
     *    AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
C
C     PASSED
C
C      PRINT *,' AFTER NEW HTOPER, ',GHAT(LMAD(1,1),1,1),
C    *   GHAT(LMAD(1,1),1,2),GHAT(LMAD(1,1),1,3),GHAT(LMAD(1,1),1,4),
C    *        GHATP(LMAD(1,1))
C      WRITE(6,*)'AT 11',GHAT(IMAX,1,1)
      IF(JITER.GT.1)THEN
       DO I=1,NCEF4S
        GHAT(I,1,1)=GHAT(I,1,1)+XHAT(I,1,1)
       ENDDO
C-CRA                 GHATP=GHATP+XHATP
C       REAL GHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          GHATP(ITMP)=GHATP(ITMP)+XHATP(ITMP)
          ENDDO
      ENDIF
C-CRA                XHAT=0.
C       REAL XHAT((JCAP+1)*(JCAP+2),NSIG,4)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG*4
          XHAT(ITMP,1,1)=0.
          ENDDO
C-CRA                XHATP=0.
C       REAL XHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          XHATP(ITMP)=0.
          ENDDO
      DO 956 I=1,NCEF4S
       PHAT(I,1,1)=-GHAT(I,1,1)
 956  CONTINUE
C-CRA                PHATP=-GHATP
C       REAL PHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          PHATP(ITMP)=-GHATP(ITMP)
          ENDDO
C        FHAT USED NOW AS WORK ARRAY
C----------
C--------
C-------- BEGINNING OF ITERATION
C--------
      RRMOLD=SDOT(4*NCEFS,GHAT(1,1,1),1,GHAT(1,1,1),1)
     *         +SDOT(NCEF,GHATP,1,GHATP,1)
      IF(JITER.EQ.1) THEN
       RRM0=RRMOLD
      END IF
      WRITE(6,2677)RRMOLD
2677  FORMAT(' DYNAMICS+MOISTURE RORIGINAL=',E12.5)
      A=.5
C---------------TEST FOR NO DATA
      PRINT *,'BEFORE TEST FOR NO DATA'
      IF(RRMOLD.LT.1.E-10) GO TO 3070
      DO 3000 ITER=1,NITER
C----------                         1/2 T -1   1/2
C---------- APPLY OPERATOR A = I + B   H O  H B
C----------
        DO I=1,NCEF4S
         FHAT(I,1,1)=PHAT(I,1,1)
        END DO
C-CRA                  FHATP=PHATP
C       REAL FHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          FHATP(ITMP)=PHATP(ITMP)
          ENDDO
C      PRINT *,' BEFORE NEW HOPER, ',FHAT(LMAD(1,1),1,1),
C    *   FHAT(LMAD(1,1),1,2),FHAT(LMAD(1,1),1,3),FHAT(LMAD(1,1),1,4),
C    *        FHATP(LMAD(1,1))
        CALL HOPER(FHAT(1,1,1),FHAT(1,1,2),FHAT(1,1,3),
     *    FHAT(1,1,4),FHATP,RU,RV,RVORT,RT,RP,RPLON,RPLAT,RQ,
     *    BHALF,BHALFP,NSIG,JCAP,NLON,NLATH,DEL2,
     *    PLN,QLN,RLN,TRIGS,IFAX,
     *    AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
C
C     PASSED
C
C     PRINT *,' AFTER NEW HOPER, ',RU(1),RV(1),RVORT(1),
C    *   RT(1),RP(1),RPLON(1),RPLAT(1),RQ(1)
C-------------------------------------ADD DIV-TEND PENALTY HERE
       IF(AMPDIVT.GT.0.) THEN
C     PRINT *,' BEFORE NEW QOPER, ',RU(1),RV(1),RVORT(1),
C    *   RT(1),RPLON(1),RPLAT(1)
        CALL QOPER(RU,RV,RVORT,RT,RPLON,RPLAT,
     *    NSIG,JCAP,NLON,NLATH,PLN,QLN,RLN,TRIGS,IFAX,
     *    DEL2,WGTS,A3,SIGL,SIGI,FHAT,ISCRA3,RLATS,
     *    RUS,RVS,RTS,RVORTS,RPLONS,RPLATS)
C      PRINT *,' AFTER NEW QOPER, ',FHAT(LMAD(1,1),1,1)
        DO K=1,NSIG
         DO I=1,NCEF
          FHAT(I,K,1)=BDIVT(IN(I),K)*FHAT(I,K,1)
         END DO
C-------------------MULTIPLY ZONAL TERMS BY EXTRA FACTOR OF 2
C------------- FOR HOMOGENEOUS, ISOTROPIC COVARIANCE IN GRID SPACE.
C-CRA                   FHAT(1:NCEF,K,1)=FACTORI(1:NCEF)*FHAT(1:NCEF,K,1)
          DO ITMP=1,NCEF
                   FHAT(ITMP,K,1)=FACTORI(ITMP)*FHAT(ITMP,K,1)
          END DO
        END DO
       END IF
C       PEND=SDOT(3*NCEFS,XHAT(1,1,1),1,XHAT(1,1,1),1)
C       PENQ=SDOT(NCEFS,XHAT(1,1,4),1,XHAT(1,1,4),1)
C       WRITE(72,*)PEND,PENQ
C       PRINT *,' INITIAL ERROR PENALTIES ',PEND,PENQ
C-CRA                  ST=RT
C       REAL ST((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          ST(ITMP)=RT(ITMP)
          ENDDO
        CALL INTPS(RP,NLATH,NLON,NPRECS,PSFILE)
        CALL INTW(RU,RV,NLATH,NLON,NSIG,NWRECS,UVFILE)
        CALL INTT(ST,NTRECS,NLATH,NLON,NSIG,TFILE)
        CALL INTQPW(RQ,NLATH,NLON,NSIG,NQRECS,NPWRECS,
     *     PWCON,QFILE,PWFILE)
C
C       PASSED
C
        CALL SATOP4(MSAT,NLATH,NLON,NSIG,RT,CSHAT,
     *    PLN,TRIGS,IFAX,JCAP,ISATV,SFILE)
C       PASSED
C-CRA                  RT=RT+ST
C       REAL RT((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RT(ITMP)=RT(ITMP)+ST(ITMP)
          ENDDO
        IF(AMPDIVT.GT.0.) THEN
         PRINT *,' BEFORE NEW QTOPER, ',FHAT(LMAD(1,1),1,1)
         CALL QTOPER(RU,RV,RVORT,RT,RPLON,RPLAT,
     *     NSIG,JCAP,NLON,NLATH,PLN,QLN,RLN,TRIGS,IFAX,
     *     DEL2,WGTS,A3,SIGL,SIGI,FHAT,ISCRA3,RLATS,
     *     RUS,RVS,RTS,RVORTS,RPLONS,RPLATS)
C
C
         PRINT *,' AFTER NEW QTOPER, ',RU(1),RV(1),RVORT(1),
     *      RT(1),RPLON(1),RPLAT(1)
        ELSE
C-CRA                   RVORT=0.
C          REAL RVORT((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RVORT(ITMP)=0.
          ENDDO
C-CRA                   RPLON=0.
C          REAL RPLON((2*NLATH+1)*(NLON+2))
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          RPLON(ITMP)=0.
          ENDDO
C-CRA                   RPLAT=0.
C          REAL RPLAT((2*NLATH+1)*(NLON+2))
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          RPLAT(ITMP)=0.
          ENDDO
        END IF
C
C
      PRINT *,' BEFORE NEW HTOPER, ',RU(1),RV(1),RVORT(1),
     *   RT(1),RP(1),RPLON(1),RPLAT(1),RQ(1)
        CALL HTOPER(FHAT(1,1,1),FHAT(1,1,2),FHAT(1,1,3),
     *    FHAT(1,1,4),FHATP,RU,RV,RVORT,RT,RP,RPLON,RPLAT,RQ,
     *    BHALF,BHALFP,NSIG,JCAP,NLON,NLATH,DEL2,
     *    PLN,QLN,RLN,TRIGS,IFAX,
     *    AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
       PRINT *,' AFTER NEW HTOPER, ',FHAT(LMAD(1,1),1,1),
     *   FHAT(LMAD(1,1),1,2),FHAT(LMAD(1,1),1,3),FHAT(LMAD(1,1),1,4),
     *        FHATP(LMAD(1,1))
C      WRITE(6,*)'AT 21',FHAT(IMAX,1,1)
        DO 958 I=1,NCEF4S
 958    FHAT(I,1,1)=PHAT(I,1,1)+FHAT(I,1,1)
C-CRA                  FHATP=PHATP+FHATP
C       REAL FHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          FHATP(ITMP)=PHATP(ITMP)+FHATP(ITMP)
          ENDDO
C      WRITE(6,*)'AT 22',FHAT(IMAX,1,1)
        A=0.
        PTGP=SDOT(4*NCEFS,FHAT(1,1,1),1,PHAT(1,1,1),1)
     *      +SDOT(NCEF,FHATP,1,PHATP,1)
        IF(PTGP .GT. 1.E-16)A=RRMOLD/PTGP
        DO 858 I=1,4*NCEFS
        XHAT(I,1,1)=XHAT(I,1,1)+A*PHAT(I,1,1)
 858    GHAT(I,1,1)=GHAT(I,1,1)+A*FHAT(I,1,1)
C-CRA                  XHATP=XHATP+A*PHATP
C       REAL XHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          XHATP(ITMP)=XHATP(ITMP)+A*PHATP(ITMP)
          ENDDO
C-CRA                  GHATP=GHATP+A*FHATP
C       REAL GHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          GHATP(ITMP)=GHATP(ITMP)+A*FHATP(ITMP)
          ENDDO
        RRMNEW=SDOT(4*NCEFS,GHAT(1,1,1),1,GHAT(1,1,1),1)
     *           +SDOT(NCEF,GHATP,1,GHATP,1)
        B=0.
        IF(ABS(RRMOLD).GT.1.E-16) B=RRMNEW/RRMOLD
        WRITE(6,2678)ITER,RRMNEW,A,B
2678    FORMAT(' DYNAMICS ITER,RNEW,A,B=',I3,3E12.5)
        IF(ITER .EQ. NITER)GO TO 3100
        IF(RRMNEW.LT.1.E-7*RRM0)GO TO 3050
        DO 608 I=1,NCEFS
          PHAT(I,1,1)=-GHAT(I,1,1)+B*PHAT(I,1,1)
          PHAT(I,1,2)=-GHAT(I,1,2)+B*PHAT(I,1,2)
          PHAT(I,1,3)=-GHAT(I,1,3)+B*PHAT(I,1,3)
          PHAT(I,1,4)=-GHAT(I,1,4)+B*PHAT(I,1,4)
608     CONTINUE
C-CRA                  PHATP=-GHATP+B*PHATP
C       REAL PHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          PHATP(ITMP)=-GHATP(ITMP)+B*PHATP(ITMP)
          ENDDO
        RRMOLD=RRMNEW
      PRINT *,'ITERATION ',ITER,' COMPLETED'
3000  CONTINUE
      GO TO 3100
3050  CONTINUE
      WRITE(6,3060)
3060  FORMAT(' ITERATION STOPPED BECAUSE RESIDUAL REDUCED BY ',
     *    'MORE THAN 7 ORDERS OF MAGNITUDE.')
      GO TO 3100
3070  CONTINUE
       WRITE(6,3075)
3075   FORMAT('   APPARENTLY NO DATA')
3100  CONTINUE
C
C--------
C-------- CONVERT TO MODEL VARIABLES
C--------
C-CRA                  GHAT=XHAT
C       REAL GHAT((JCAP+1)*(JCAP+2),NSIG,4)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG*4
          GHAT(ITMP,1,1)=XHAT(ITMP,1,1)
          ENDDO
C-CRA                  GHATP=XHATP
C       REAL GHATP((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          GHATP(ITMP)=XHATP(ITMP)
          ENDDO
        CALL HOPERS(GHAT(1,1,1),GHAT(1,1,2),GHAT(1,1,3),
     *    GHAT(1,1,4),GHATP,BHALF,BHALFP,NSIG,JCAP,
     *     AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
C-CRA                  FHAT=GHAT
C       REAL FHAT((JCAP+1)*(JCAP+2),NSIG,4)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG*4
          FHAT(ITMP,1,1)=GHAT(ITMP,1,1)
          ENDDO
C-CRA                  PS=GHATP
C       REAL PS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          PS(ITMP)=GHATP(ITMP)
          ENDDO
C--------
C-------- 14. ADD INCREMENT TO GUESS AND WRITE OUT
C--------
C--------
C-------- READ IN GUESS, PUTTING INTO INTERNAL FORMAT.
C-------- SCRA AND AQR USED AS SCRATCH
C--------
      CALL RDGESC(PHAT(1,1,1),PHAT(1,1,2),PHAT(1,1,3),
     *  PHAT(1,1,4),FACTOR,FACTORI,HOURG,
     *  IDATEG,SIGI,SIGL,INEXT,JCAP,NSIG,ON85,
     *  ML2LM,FACTSLM,FACTVLM)
C--------
C-------- ADD INCREMENTS
C--------
      DO L=1,NSIG*(JCAP+1)*(JCAP+2)
       FHAT(L,1,1)=FHAT(L,1,1)+PHAT(L,1,1)
       FHAT(L,1,2)=FHAT(L,1,2)+PHAT(L,1,2)
       FHAT(L,1,3)=FHAT(L,1,3)+PHAT(L,1,3)
       FHAT(L,1,4)=FHAT(L,1,4)+PHAT(L,1,4)
      END DO
C-CRA                PS=FACTOR+PS
C       REAL PS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          PS(ITMP)=FACTOR(ITMP)+PS(ITMP)
          ENDDO
C--------
C-------- GET GUESS Q TO USE AS REFERENCE WHEN LIMITING ANALYSIS Q
C--------   (ONLY AFTER LAST OUTER ITERATION)
C-------
      IF (JITER.EQ.MITER) THEN
       CALL RDGESC(PHAT(1,1,1),PHAT(1,1,2),PHAT(1,1,3),
     *   PHAT(1,1,4),FACTOR,FACTORI,HOURG,
     *   IDATEG,SIGI,SIGL,INGES,JCAP,NSIG,ON85,
     *   ML2LM,FACTSLM,FACTVLM)
       CALL S2G0(PS,RV,JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
       DO K=1,NSIG
        CALL S2G0(PHAT(1,K,4),RU((K-1)*NGRP+1),JCAP,NLON,
     *     NLATH,PLN,TRIGS,IFAX)
        CALL S2G0(FHAT(1,K,3),RT((K-1)*NGRP+1),JCAP,NLON,
     *     NLATH,PLN,TRIGS,IFAX)
       END DO
C-CRA                 RQ=RU
C       REAL RQ((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RQ(ITMP)=RU(ITMP)
          ENDDO
       CALL GENQSAT(RT,RQ,NLATH,NLON,NSIG,
     *    RV,SIGL)
C-CRA                 RV=0.
C       REAL RV((2*NLATH+1)*(NLON+2)*NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RV(ITMP)=0.
          ENDDO
       DO I=1,NGRPS
        IF(RQ(I) .GT. 0. .AND. RU(I) .GT. RQ(I))
     *      RQ(I)=RU(I)
       END DO
       DO I=1,NGRPS                         
        IF(RU(I) .LT. 0.  )RV(I)=RU(I)
       END DO
       DO K=1,NSIG
        CALL S2G0(FHAT(1,K,4),RU((K-1)*NGRP+1),JCAP,NLON,
     *     NLATH,PLN,TRIGS,IFAX)
       END DO
       NSUPER=0
       NNQ=0
       DO I=1,NGRPS
         IF(RU(I) .GT. 0. .AND. RU(I) .GT. RQ(I)) THEN
          RU(I)=RQ(I)
          NSUPER=NSUPER+1
         END IF
C       END IF
       END DO
       DO I=1,NGRPS
C       IF(RU(I) .LT. 0. .AND. RU(I) .LT. RQ(I)) THEN
C        RU(I)=0.
        IF(RU(I) .LT. 0. .AND. RU(I) .LT. RV(I)) THEN
         RU(I)=MIN(0.,RV(I))
         NNQ=NNQ+1
        END IF
       END DO
       PRINT *,' NUMBER OF SUPERSATURATED POINTS = ',NSUPER
       PRINT *,' NUMBER OF NEGATIVE Q POINTS = ',NNQ
       DO K=1,NSIG
        CALL G2S0(FHAT(1,K,4),RU((K-1)*NGRP+1),JCAP,NLON,NLATH,
     *           WGTS,PLN,TRIGS,IFAX)
       END DO
      END IF
      HOURG=0.
C--------
C-------- WRITE OUT ANALYSIS
C--------
      CALL WRANLC(FHAT(1,1,1),FHAT(1,1,2),FHAT(1,1,3),
     *  FHAT(1,1,4),
     *  PS,FACTORI,HOURG,IDATEG,SIGI,SIGL,IOANL,
     *   JCAP,NSIG,ON85,ON85DT,LM2ML,FACTSML,FACTVML)
      RETURN
      END
