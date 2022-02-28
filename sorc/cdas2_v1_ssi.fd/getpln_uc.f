       SUBROUTINE GETPLN(PLN,QLN,RLN,JCAP,NLATH,AP,BP,SLAT,PE0,
     *          QE0,RO0,AQR,BQR,GR,CLAT,DEL2,DEL2OUT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    GETPLN     GENERATE LEGENDRE POLYNOMIALS
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-09-21
C
C ABSTRACT: SUMMATION OF SCALAR SPHERICAL HARMONIC SERIES.
C
C PROGRAM HISTORY LOG:
C   90-09-21  PARRISH
C
C   INPUT ARGUMENT LIST:
C     JCAP     - TRIANGULAR TRUNCATION
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     AP,BP    - RECURSION CONSTANTS FOR SPHERICAL HARMONICS
C     SLAT     - SIN(GAUSSIAN LATITUDES)
C     PE0      - STARTING FUNCTIONS FOR SPHERICAL HARMONICS
C
C   OUTPUT ARGUMENT LIST:
C     PLN      - LEGENDRE POLYNOMIALS
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA             DIMENSION AP(0:JCAP,0:JCAP),AQR(0:JCAP,0:JCAP)
C-CRA             DIMENSION BP(0:JCAP,0:JCAP),BQR(0:JCAP,0:JCAP)
C-CRA             DIMENSION GR(0:JCAP,0:JCAP)
C-CRA             DIMENSION DEL2(0:JCAP,0:JCAP)
C-CRA             DIMENSION SLAT(NLATH),CLAT(NLATH)
C-CRA             DIMENSION PE0(NLATH,0:JCAP)
C-CRA             DIMENSION QE0(NLATH,0:JCAP)
C-CRA             DIMENSION RO0(NLATH,0:JCAP)
C-CRA             DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION DEL2OUT((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION IADR(0:JCAP,0:JCAP)
C-CRA             REAL PE(NLATH,0:JCAP),PO(NLATH,0:JCAP)
C-CRA             REAL QE(NLATH,0:JCAP),QO(NLATH,0:JCAP)
C-CRA             REAL RE(NLATH,0:JCAP),RO(NLATH,0:JCAP)
 
             DIMENSION AP(0:62,0:62),AQR(0:62,0:62)
             DIMENSION BP(0:62,0:62),BQR(0:62,0:62)
             DIMENSION GR(0:62,0:62)
             DIMENSION DEL2(0:62,0:62)
             DIMENSION SLAT(48),CLAT(48)
             DIMENSION PE0(48,0:62)
             DIMENSION QE0(48,0:62)
             DIMENSION RO0(48,0:62)
             DIMENSION PLN((62+1)*(62+2),48)
             DIMENSION QLN((62+1)*(62+2),48)
             DIMENSION RLN((62+1)*(62+2),48)
             DIMENSION DEL2OUT((62+1)*(62+2))
             DIMENSION IADR(0:62,0:62)
             REAL PE(48,0:62),PO(48,0:62)
             REAL QE(48,0:62),QO(48,0:62)
             REAL RE(48,0:62),RO(48,0:62)
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C--------
         II=-1
         DO M=0,JCAP
          DO L=0,JCAP-M
           II=II+2
           IADR(L,M)=II
          END DO
         END DO
         DO M=0,JCAP
          DEL2OUT(IADR(0,M))=DEL2(M,0)
          DEL2OUT(IADR(0,M)+1)=0.
          IF(M.LT.JCAP) THEN
           DO L=1,JCAP-M
            DEL2OUT(IADR(L,M))=DEL2(M,L)
            DEL2OUT(IADR(L,M)+1)=DEL2(M,L)
           END DO
          END IF
         END DO
C-CRA                   PLN=0.
C          DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NLATH
          PLN(ITMP,1)=0.
          ENDDO
C-CRA                   QLN=0.
C          DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NLATH
          QLN(ITMP,1)=0.
          ENDDO
C-CRA                   RLN=0.
C          DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NLATH
          RLN(ITMP,1)=0.
          ENDDO
C-CRA                   PO=0.
C          REAL PE(NLATH,0:JCAP),PO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          PO(ITMP,JTMP)=0.
          ENDDO
          ENDDO
C-CRA                   PE=PE0
C          REAL PE(NLATH,0:JCAP),PO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          PE(ITMP,JTMP)=PE0(ITMP,JTMP)
          ENDDO
          ENDDO
C-CRA                   QO=0.
C          REAL QE(NLATH,0:JCAP),QO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          QO(ITMP,JTMP)=0.
          ENDDO
          ENDDO
C-CRA                   QE=QE0
C          REAL QE(NLATH,0:JCAP),QO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          QE(ITMP,JTMP)=QE0(ITMP,JTMP)
          ENDDO
          ENDDO
C-CRA                   RE=0.
C          REAL RE(NLATH,0:JCAP),RO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          RE(ITMP,JTMP)=0.
          ENDDO
          ENDDO
C-CRA                   RO=RO0
C          REAL RE(NLATH,0:JCAP),RO(NLATH,0:JCAP)
          DO JTMP=0,JCAP
          DO ITMP=1,NLATH
          RO(ITMP,JTMP)=RO0(ITMP,JTMP)
          ENDDO
          ENDDO
         DO L=0,JCAP
          DO M=0,JCAP-L,2
C------------ FIRST EVEN TERMS (M=0,2,...)
           DO J=1,NLATH
            PLN(IADR(L,M),J)=PE(J,L)
            QLN(IADR(L,M),J)=QE(J,L)
            RLN(IADR(L,M),J)=RO(J,L)
           END DO
C------------ NOW DO ODD  (M=1,3,...)
           IF(M+1.LE.JCAP-L) THEN
            MP=M+1
            DO J=1,NLATH
              PO(J,L)=AP(M,L)*SLAT(J)*PE(J,L)+BP(M,L)*PO(J,L)
              QO(J,L)=AQR(M,L)*SLAT(J)*QE(J,L)
     *                     +BQR(M,L)*QO(J,L)
              RE(J,L)=AQR(M,L)*SLAT(J)*RO(J,L)
     *                 +BQR(M,L)*RE(J,L)+GR(M,L)*PE(J,L)*CLAT(J)
            END DO
            DO J=1,NLATH
             PLN(IADR(L,MP),J)=PO(J,L)
             QLN(IADR(L,MP),J)=QO(J,L)
             RLN(IADR(L,MP),J)=RE(J,L)
            END DO
C-------------- GET NEXT PE
            DO J=1,NLATH
             PE(J,L)=AP(MP,L)*SLAT(J)*PO(J,L)+BP(MP,L)*PE(J,L)
             QE(J,L)=AQR(MP,L)*SLAT(J)*QO(J,L)
     *                 +BQR(MP,L)*QE(J,L)
             RO(J,L)=AQR(MP,L)*SLAT(J)*RE(J,L)
     *                 +BQR(MP,L)*RO(J,L)+GR(MP,L)*PO(J,L)*CLAT(J)
            END DO
           END IF
          END DO
         END DO
         DO I=1,(JCAP+1)*(JCAP+2)*NLATH,2
          PLN(I+1,1)=PLN(I,1)
          QLN(I+1,1)=QLN(I,1)
          RLN(I+1,1)=RLN(I,1)
         END DO
       RETURN
       END
