      SUBROUTINE GENQSAT(T,QSAT,NLATH,NLON,NSIG,PS,SLG)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    GENQSAT     OBTAIN QSAT FOR GIVEN T.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-11
C
C ABSTRACT: OBTAIN SATURATION SPECIFIC HUMIDITY FOR GIVEN TEMPERATURE.
C
C PROGRAM HISTORY LOG:
C   90-10-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     T,PS     - TEMPERATURE, LOG(PSFC) ON GAUSSIAN GRID
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     NLON     - NUMBER OF LONGITUDES
C     NSIG     - NUMBER OF SIGMA LEVELS
C     SLG      - SIGMA VALUES AT MID-POINT OF MODEL LEVELS
C
C   OUTPUT ARGUMENT LIST:
C     QSAT     - SATURATION SPECIFIC HUMIDITY 
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
      DIMENSION T(2*NLATH+1,NLON+2,NSIG)
      DIMENSION QSAT(2*NLATH+1,NLON+2,NSIG)
      DIMENSION PS(2*NLATH+1,NLON+2),SLG(NSIG)
      REAL L0
      EPS=.622
      CP=1005.
      CL=4187.
      CPV=1876.5
      RV=461.5
      L0=2.501E6
      T0=273.16
      ES0=611.0
      FACT1 = (CPV - CL) / RV
      FACT2 = (L0 + (CL - CPV) * T0) / RV
      FACT3 = 1. / T0
      OMEPS=1.-EPS
C--------
      DO 200 K=1,NSIG
        DO 100 I=1,NLON
          DO 100 J=1,2*NLATH
            PW=EXP(PS(J,I))*SLG(K)
            IF(QSAT(J,I,K) .LT. 0.)QSAT(J,I,K)=0.
            TV=T(J,I,K)/(1.+0.61*QSAT(J,I,K))
            IF(PW.LT.5.) PW=5.
            PW=1000.*PW
            ES = ES0 * (TV / T0) ** FACT1 *
     1          EXP ( FACT2 * (FACT3 - 1. / TV))
            QS = EPS * ES / (PW - OMEPS * ES)
            IF(QS .LT. QSAT(J,I,K))THEN
            TV=T(J,I,K)/(1.+0.61*QS)
            ES = ES0 * (TV / T0) ** FACT1 *
     1          EXP ( FACT2 * (FACT3 - 1. / TV))
            QS = EPS * ES / (PW - OMEPS * ES)
            TV=T(J,I,K)/(1.+0.61*QS)
            ES = ES0 * (TV / T0) ** FACT1 *
     1          EXP ( FACT2 * (FACT3 - 1. / TV))
            QS = EPS * ES / (PW - OMEPS * ES)
            END IF
            QSAT(J,I,K)=QS
100       CONTINUE
200   CONTINUE
      RETURN
      END
