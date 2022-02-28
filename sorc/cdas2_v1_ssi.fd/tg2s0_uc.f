       SUBROUTINE TG2S0(TS,T,JCAP,NLON,NLATH,WGTS,PLN,TRIGS,IFAX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    TRANSPOSE OF G2S0
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-09-21
C
C ABSTRACT: SUMMATION OF SCALAR SPHERICAL HARMONIC SERIES.
C
C PROGRAM HISTORY LOG:
C   90-09-21  PARRISH
C
C   INPUT ARGUMENT LIST:
C     TS       - SPECTRAL COEFS
C     JCAP     - TRIANGULAR TRUNCATION
C     NLON     - NUMBER OF LONGITUDES
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     PLN      - SPHERICAL HARMONICS
C     TRIGS,IFAX - USED BY FFT
C
C   OUTPUT ARGUMENT LIST:
C     T        - VALUES OF DESIRED FIELD ON GAUSSIAN GRID
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA             DIMENSION TS((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION T(2*NLATH+1,NLON+2)
C-CRA             DIMENSION TRIGS(NLON*2),IFAX(10)
C-CRA             DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION WGTS(NLATH)
C-CRA             DIMENSION WORK(2*(2*NLATH+1)*(NLON+2))
C-CRA             DIMENSION TE(2*JCAP+2),TO(2*JCAP+2)
C-CRA             DIMENSION FACTOR(2*JCAP+2,NLATH)
 
             DIMENSION TS((62+1)*(62+2))
             DIMENSION T(2*48+1,192+2)
             DIMENSION TRIGS(192*2),IFAX(10)
             DIMENSION PLN((62+1)*(62+2),48)
             DIMENSION WGTS(48)
             DIMENSION WORK(2*(2*48+1)*(192+2))
             DIMENSION TE(2*62+2),TO(2*62+2)
             DIMENSION FACTOR(2*62+2,48)
C--------
         DO J=1,NLATH
          FACTOR(1,J)=WGTS(J)/NLON
          FACTOR(2,J)=0.
          DO I=3,2*JCAP+2
           FACTOR(I,J)=.5*FACTOR(1,J)
          END DO
         END DO
C-CRA                   T=0.
C          DIMENSION T(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          T(ITMP,1)=0.
          ENDDO
         DO J=1,NLATH
          JR=2*NLATH+1-J
          II0=0
C-CRA                    TE=0.
C          DIMENSION TE(2*JCAP+2),TO(2*JCAP+2)
          DO ITMP=1,2*JCAP+2
          TE(ITMP)=0.
          ENDDO
C-CRA                    TO=0.
C          DIMENSION TE(2*JCAP+2),TO(2*JCAP+2)
          DO ITMP=1,2*JCAP+2
          TO(ITMP)=0.
          ENDDO
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            TE(LL)=TE(LL)+PLN(II0+LL,J)*TS(II0+LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             TO(LL)=TO(LL)+PLN(II0+LL,J)*TS(II0+LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
C----------
C---------- NOW COMBINE EVEN AND ODD PARTS
C----------
          DO LL=1,2*(JCAP+1)
           T(J,LL)=(TE(LL)+TO(LL))*FACTOR(LL,J)
           T(JR,LL)=(TE(LL)-TO(LL))*FACTOR(LL,J)
          END DO
         END DO
C--------
C-------- FINALLY DO FOURIER SUMS IN LONGITUDE
C--------
         LOT=NLATH*2
         NLAX=LOT+1
C-CRA             CALL RFFTMLT(T,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
             CALL FFT99M (T,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
       RETURN
       END
