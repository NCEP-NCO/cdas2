       SUBROUTINE S2MG2X(TS,T,JCAP,NLATH,NLON,PLN)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    S2MG2X     SPECTRAL VARIANCE TO GRID VARIANCE
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-12-07
C
C ABSTRACT: COMPUTE GRID VARIANCE FROM SPECTRAL DIAG COVARIANCE.
C PROGRAM HISTORY LOG:
C   90-12-07  PARRISH
C
C   INPUT ARGUMENT LIST:
C     TS       - SPECTRAL VARIANCE SPECTRUM
C     JCAP     - TRIANGULAR TRUNCATION
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     NLON     - NUMBER OF LONGITUDES
C     SLAT     - SIN(GAUSSIAN LATITUDES)
C     PLN      - SPHERICAL HARMONICS
C
C   OUTPUT ARGUMENT LIST:
C     T        - VARIANCE ON GRID
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA             DIMENSION TS((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION T(2*NLATH+1,NLON+2)
C-CRA             DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION WORK(2*NLATH+1,NLON+2)
C-CRA             DIMENSION COS2(0:JCAP,NLON),SIN2(0:JCAP,NLON)
C-CRA             DIMENSION TE(2*JCAP+2),TO(2*JCAP+2)
 
             DIMENSION TS((62+1)*(62+2))
             DIMENSION T(2*48+1,192+2)
             DIMENSION PLN((62+1)*(62+2),48)
             DIMENSION WORK(2*48+1,192+2)
             DIMENSION COS2(0:62,192),SIN2(0:62,192)
             DIMENSION TE(2*62+2),TO(2*62+2)
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C--------
         DLON=8.*ATAN(1.)/NLON
         DO I=1,NLON
          ANGLE=(I-1.)*DLON
          COS2(0,I)=1.
          SIN2(0,I)=0.
          DO L=1,JCAP
           ARG=L*ANGLE
           COS2(L,I)=4.*COS(ARG)**2
           SIN2(L,I)=4.*SIN(ARG)**2
          END DO
         END DO
C--------
C-------- NOW SUM IN LATITUDE (AFTER ZEROING OUTPUT ARRAYS)
C--------
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
            TE(LL)=TE(LL)+PLN(II0+LL,J)**2*TS(II0+LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             TO(LL)=TO(LL)+PLN(II0+LL,J)**2*TS(II0+LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
C----------
C---------- NOW COMBINE EVEN AND ODD PARTS
C----------
          DO LL=1,2*(JCAP+1)
           T(J,LL)=TE(LL)+TO(LL)
           T(JR,LL)=TE(LL)+TO(LL)
          END DO
         END DO
C--------
C------- FINALLY DO SQUARED FOURIER SUMS IN LONGITUDE
C-------
C-CRA                   WORK=T
C          DIMENSION WORK(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          WORK(ITMP,1)=T(ITMP,1)
          ENDDO
C-CRA                   T=0.
C          DIMENSION T(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          T(ITMP,1)=0.
          ENDDO
         DO L=0,JCAP
          LR=2*L+1
          LI=2*L+2
          DO I=1,NLON
           DO J=1,2*NLATH
            T(J,I)=T(J,I)+WORK(J,LR)*COS2(L,I)+WORK(J,LI)*SIN2(L,I)
           END DO
          END DO
         END DO
       RETURN
       END
