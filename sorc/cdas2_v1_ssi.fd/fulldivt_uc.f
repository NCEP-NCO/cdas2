      SUBROUTINE FULLDIVT(U,V,T,VORT,DIV,PLON,PLAT,Z0,
     *       NSIG,JCAP,NLON,NLATH,PLN,QLN,RLN,TRIGS,IFAX,DEL2,
     *       WGTS,A3,SIGL,SIGI,JITER,DS,RLATS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FULLDIVT    COMPUTE FULL FIELD DIVERGENCE TENDENCY
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 94-02-11
C
C ABSTRACT: GET DIVERGENCE TENDENCY (ONLY DYNAMIC TERMS SO FAR)
C
C PROGRAM HISTORY LOG:
C   94-02-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     U,V,T    - U, V, T ON GAUSSIAN GRID
C     VORT,DIV - VORT, DIV ON GAUSSIAN GRID
C     PLON,PLAT- LON, LAT DERIVATIVES OF PSFC ON GAUSSIAN GRID
C     Z0       - TERRAIN HEIGHT ON GAUSSIAN GRID
C     NSIG     - NUMBER OF SIGMA LEVELS
C     JCAP     - TRIANGULAR TRUNCATION OF SPECTRAL COEFFICIENTS
C     NLON     - NUM OF LONGITUDES
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     RLATS    - GAUSSIAN LATITUDES
C     TRIGS,IFAX  - USED BY FFT
C     DEL2     - CONSTANTS FOR APPLICATION OF DEL**2 OPERATOR
C     WGTS     - GAUSSIAN INTEGRATION WEIGHTS.
C     A3       - HYDROSTATIC MATRIX
C     SIGL,SIGI- SIGMA COORDINATES
C     JITER    - OUTER ITERATION COUNTER
C
C   OUTPUT ARGUMENT LIST:
C     DIV      - ON OUTPUT, CONTAINS SIGDOT, THE VERTICAL VELOCITY
C     DS       - DIV TENDENCY (SPECTRAL COEFFICIENTS)
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C-----------
C-CRA          DIMENSION SIGL(NSIG),SIGI(NSIG+1)
C-CRA          DIMENSION VORT(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION DIV(2*NLATH+1,NLON+2,NSIG)
C-CRA          REAL T(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION A3(NSIG,NSIG)
C-CRA          DIMENSION DEL2((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION WGTS(2*NLATH)
C-CRA          DIMENSION U(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION V(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION PLON(2*NLATH+1,NLON+2),PLAT(2*NLATH+1,NLON+2)
C-CRA          DIMENSION Z0(2*NLATH+1,NLON+2)
C-CRA          DIMENSION TRIGS(NLON*2),IFAX(10)
C-CRA          DIMENSION DS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION RLATS(2*NLATH)
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C-CRA          DIMENSION TS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION BIGE(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION UD(2*NLATH+1,NLON+2),VD(2*NLATH+1,NLON+2)
C-CRA          DIMENSION UW(2*NLATH+1,NLON+2,NSIG),VW(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION CORIOLIS(2*NLATH+1,NLON+2)
C-CRA          DIMENSION P(2*NLATH+1,NLON+2)
C-CRA          DIMENSION DSMS(NSIG)
C-CRA          DIMENSION FACTOR((JCAP+1)*(JCAP+2))
C-----------
          DIMENSION SIGL(28),SIGI(28+1)
          DIMENSION VORT(2*48+1,192+2,28)
          DIMENSION DIV(2*48+1,192+2,28)
          REAL T(2*48+1,192+2,28)
          DIMENSION A3(28,28)
          DIMENSION DEL2((62+1)*(62+2))
          DIMENSION WGTS(2*48)
          DIMENSION U(2*48+1,192+2,28)
          DIMENSION V(2*48+1,192+2,28)
          DIMENSION PLON(2*48+1,192+2)
          DIMENSION PLAT(2*48+1,192+2)
          DIMENSION Z0(2*48+1,192+2)
          DIMENSION TRIGS(192*2),IFAX(10)
          DIMENSION DS((62+1)*(62+2),28)
          DIMENSION PLN((62+1)*(62+2),48)
          DIMENSION QLN((62+1)*(62+2),48)
          DIMENSION RLN((62+1)*(62+2),48)
          DIMENSION RLATS(2*48)
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
          DIMENSION TS((62+1)*(62+2),28)
          DIMENSION BIGE(2*48+1,192+2,28)
          DIMENSION UD(2*48+1,192+2),VD(2*48+1,192+2)
          DIMENSION UW(2*48+1,192+2,28)
          DIMENSION VW(2*48+1,192+2,28)
          DIMENSION CORIOLIS(2*48+1,192+2)
          DIMENSION P(2*48+1,192+2)
          DIMENSION DSMS(28)
          DIMENSION FACTOR((62+1)*(62+2))
C--------
C-------- COMPUTE UD,VD, BIGE        
C--------
         II=-1
         DO M=0,JCAP
          II=II+2
          FACTOR(II)=.5
          FACTOR(II+1)=0.
          IF(M.LT.JCAP) THEN
           DO L=1,JCAP-M
            II=II+2
            FACTOR(II)=1.
            FACTOR(II+1)=1.
           END DO
          END IF
         END DO
         NG=(2*NLATH+1)*NLON
         NC=(JCAP+1)*(JCAP+2)
C-CRA             CORIOLIS=0.
             DO I=1,2*NLATH+1*NLON+2
                CORIOLIS(I,1)=0.
             ENDDO 
         OMEGA=CONMC('OMEGA$')
C        CALL EXIT
         DO J=1,NLATH
          JR=2*NLATH+1-J
C-CRA              CORIOLIS(J,1:NLON)=2.*OMEGA*SIN(RLATS(J))
C-CRA              CORIOLIS(JR,1:NLON)=-CORIOLIS(J,1:NLON)
              DO K=1,NLON
              CORIOLIS(J,K)=2.*OMEGA*SIN(RLATS(J))
              CORIOLIS(JR,K)=-CORIOLIS(J,K)
              ENDDO
         END DO
         GASCON=CONMC('RD$')
         EACCEL=9.8
C---------------------GET SIGDOT  (OVERWRITES DIVERGENCE)
C-CRA             P=0.
             DO I=1,(2*NLATH+1)*(NLON+2)
             P(I,1)=0.
             ENDDO
         DO K=1,NSIG
          DO I=1,NG
           AK=(U(I,1,K)*PLON(I,1)+V(I,1,K)*PLAT(I,1)
     *         +DIV(I,1,K))*(SIGI(K)-SIGI(K+1))
           P(I,1)=P(I,1)+AK
           DIV(I,1,K)=P(I,1)
          END DO
         END DO
         SIGSUM=0.
         DO K=1,NSIG
          SIGSUM=SIGSUM+SIGI(K)-SIGI(K+1)
          DO I=1,NG
           DIV(I,1,K)=DIV(I,1,K)-SIGSUM*P(I,1)
          END DO
         END DO
C--------------------NOW COMPUTE DIV TENDENCY (DYNAMICS ONLY)
         DO K=1,NSIG
          DO J=1,2*NLATH
           DO I=1,NLON
            TERM1=V(J,I,K)*(VORT(J,I,K)+CORIOLIS(J,I))
            TERM2=GASCON*T(J,I,K)*PLON(J,I)
            TERM3=0.
            IF(K.GT.1) 
     *        TERM3=TERM3+.5*DIV(J,I,K-1)*(U(J,I,K)-U(J,I,K-1))
     *            /(SIGL(K)-SIGL(K-1))
            IF(K.LT.NSIG)
     *        TERM3=TERM3+.5*DIV(J,I,K)*(U(J,I,K)-U(J,I,K+1))
     *            /(SIGL(K)-SIGL(K+1))
            UD(J,I)=TERM1-TERM2-TERM3
            TERM1=-U(J,I,K)*(VORT(J,I,K)+CORIOLIS(J,I))
            TERM2=GASCON*T(J,I,K)*PLAT(J,I)
            TERM3=0.
            IF(K.GT.1) 
     *        TERM3=TERM3+.5*DIV(J,I,K-1)*(V(J,I,K)-V(J,I,K-1))
     *            /(SIGL(K)-SIGL(K-1))
            IF(K.LT.NSIG)
     *        TERM3=TERM3+.5*DIV(J,I,K)*(V(J,I,K)-V(J,I,K+1))
     *            /(SIGL(K)-SIGL(K+1))
            VD(J,I)=TERM1-TERM2-TERM3
            TERM1=.5*(U(J,I,K)**2+V(J,I,K)**2)
            BIGE(J,I,K)=TERM1
           END DO
          END DO
C-CRA              UW(1:NG,1,K)=UD(1:NG,1)
C-CRA              VW(1:NG,1,K)=VD(1:NG,1)
              DO I=1,NG
              UW(I,1,K)=UD(I,1)
              VW(I,1,K)=VD(I,1)
              ENDDO
C-CRA              UD=EACCEL*Z0
              DO I=1,(2*NLATH+1)*(NLON+2)
              UD(I,1)=EACCEL*Z0(I,1)
              ENDDO
          DO L=1,NSIG
C-CRA               UD(1:NG,1)=UD(1:NG,1)+EACCEL*A3(K,L)*T(1:NG,1,L)
               DO I=1,NG
               UD(I,1)=UD(I,1)+EACCEL*A3(K,L)*T(I,1,L)
               ENDDO
          END DO
          DO J=1,2*NLATH
           DO I=1,NLON
            BIGE(J,I,K)=BIGE(J,I,K)+UD(J,I)
           END DO
          END DO
         END DO
C---------
C--------- NOW GET DIV (UD,VD) - DEL2 (BIGE)
C---------
         DO K=1,NSIG
          CALL G2S0(TS(1,K),BIGE(1,1,K),JCAP,NLON,NLATH,
     *         WGTS,PLN,TRIGS,IFAX)
          CALL GRAD2S(DS(1,K),UW(1,1,K),VW(1,1,K),JCAP,NLON,NLATH,
     *                QLN,RLN,TRIGS,IFAX,WGTS,DEL2)
C-CRA              TS(1:NC,K)=DEL2(1:NC)*TS(1:NC,K)
C-CRA              DS(1:NC,K)=DS(1:NC,K)+TS(1:NC,K)
              DO I=1,NC
              TS(I,K)=DEL2(I)*TS(I,K)
              DS(I,K)=DS(I,K)+TS(I,K)
              ENDDO
         END DO
C-CRA             DSMS=0.
             DO I=1,NSIG
             DSMS(I)=0.
             ENDDO
         DO K=1,NSIG
          DO I=1,NC
           DSMS(K)=DSMS(K)+FACTOR(I)*DS(I,K)**2
          END DO
         END DO
C        PRINT *,' FOR OUTER ITERATION = ',JITER,' DIVT STATS FOLLOW'
C        WRITE(6,50)DSMS
C50       FORMAT(1H ,6E13.4)
       RETURN
       END
