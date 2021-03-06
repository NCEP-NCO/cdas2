       SUBROUTINE TS2G0(TS,T,JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    TS2G0       TRANSPOSE OF S2G0
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
C-CRA             DIMENSION WORK(2*(2*NLATH+1)*(NLON+2))
C-CRA             DIMENSION WGTS(2*JCAP+2),TE(2*JCAP+2),TO(2*JCAP+2)
 
             DIMENSION TS((62+1)*(62+2))
             DIMENSION T(2*48+1,192+2)
             DIMENSION TRIGS(192*2),IFAX(10)
             DIMENSION PLN((62+1)*(62+2),48)
             DIMENSION WORK(2*(2*48+1)*(192+2))
             DIMENSION WGTS(2*62+2),TE(2*62+2),TO(2*62+2)
C--------
C-CRA                   WGTS=2.*NLON
C          DIMENSION WGTS(2*JCAP+2),TE(2*JCAP+2),TO(2*JCAP+2)
          DO ITMP=1,2*JCAP+2
          WGTS(ITMP)=2.*NLON
          ENDDO
         WGTS(1)=NLON
         WGTS(2)=0.
C--------
C-------- FIRST DO FOURIER ANALYSIS IN LONGITUDE
C--------
         LOT=NLATH*2
         NLAX=LOT+1
C-CRA             CALL RFFTMLT(T,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
             CALL FFT99M (T,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
C-CRA                   TS=0.
C          DIMENSION TS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          TS(ITMP)=0.
          ENDDO
         DO J=1,NLATH
          JR=2*NLATH+1-J
C---------- SEPARATE EVEN AND ODD PARTS
          DO LL=1,2*JCAP+2
           TE(LL)=(T(J,LL)+T(JR,LL))*WGTS(LL)
           TO(LL)=(T(J,LL)-T(JR,LL))*WGTS(LL)
          END DO
          II0=0
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            TS(II0+LL)=TS(II0+LL)+PLN(II0+LL,J)*TE(LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             TS(II0+LL)=TS(II0+LL)+PLN(II0+LL,J)*TO(LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
         END DO
       RETURN
       END
