       SUBROUTINE SATC(T,NSIG,JCAP,NLON,NLATH,
     *            PLN,TRIGS,IFAX,CSHAT,ISATV)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    SATC      SAT. ERROR CORRELATION OPERATOR
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 91-08-15
C
C ABSTRACT: OPERATES SATELLITE ERROR CORRELATION ON FIELDS
C
C PROGRAM HISTORY LOG:
C   91-08-15  PARRISH
C
C   INPUT ARGUMENT LIST:
C     T        - GRIDDED INPUT TEMPERATURE FIELD
C     NSIG     - NUMBER OF SIGMA LEVELS
C     JCAP     - TRIANGULAR TRUNCATION
C     NLON     - NUMBER OF LONGITUDES
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     PLN      - SPHERICAL HARMONICS
C     TRIGS,IFAX - USED BY FFT
C     CSHAT    - SPECTRAL WEIGHTS
C     ISATV    - FLAGS FOR INDICATING WHETHER OR NOT TO APPLY 
C                CORRELATION, >0 YES, <0 NO
C
C   OUTPUT ARGUMENT LIST:
C     T        - GRIDDED OUTPUT TEMPERATURE FIELD
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION CSHAT((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION T(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION TRIGS(NLON*2),IFAX(10),ISATV(NSIG)
C-CRA          DIMENSION TS((JCAP+1)*(JCAP+2))
 
          DIMENSION CSHAT((62+1)*(62+2))
          DIMENSION PLN((62+1)*(62+2),48)
          DIMENSION T(2*48+1,192+2,28)
          DIMENSION TRIGS(192*2),IFAX(10),ISATV(28)
          DIMENSION TS((62+1)*(62+2))
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C--------
      NBSIG=1
      NESIG=NSIG
      DO K=2,NSIG
      IF(ISATV(K-1) .LT. 0.)THEN
      IF(NBSIG .EQ. K-1)NBSIG=K
      END IF
      END DO
      DO K=NSIG-1,1,-1
      IF(ISATV(K+1) .LT. 0.)THEN
      IF(NESIG .EQ. K+1)NESIG=K
      END IF
      END DO
      NSIGL=NESIG-NBSIG+1
      IF(NSIGL .LE. 0)THEN
      PRINT *,' NEGATIVE SIGMA LEVELS IN SATC'
      STOP
      END IF
      DO K=NBSIG,NESIG
       CALL TS2G0(TS,T(1,1,K),JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
C-----------------------------MULTIPLY COEFS BY CORRELATION
          DO I=1,(JCAP+1)*(JCAP+2)
           TS(I)=TS(I)*CSHAT(I)
          END DO
       CALL S2G0(TS,T(1,1,K),JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
      END DO
      RETURN
      END