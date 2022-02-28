      SUBROUTINE M1RCONS(AP,BP,AQR,BQR,GR,DEL2,JCAP)
C...TRANSLATED BY FPP 6.0 (3.06G3) 08/09/95  14:56:24    
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    M1RCONS    COMPUTE LEGENDRE GENERATOR CONSTANTS
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-09-21
C
C ABSTRACT: GET GENERATOR CONSTANTS NEEDED FOR LEGENDRE TRANSFORMS
C
C PROGRAM HISTORY LOG:
C   90-09-21  PARRISH
C
C   INPUT ARGUMENT LIST:
C     JCAP     - TRIANGULAR TRUNCATION
C
C   OUTPUT ARGUMENT LIST:
C     AP,BP,AQR,BQR,GR - VARIOUS RECURSION CONSTANTS
C     DEL2     - N*(N+1)/(A**2)
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
      DIMENSION AP(0:JCAP,0:JCAP)
      DIMENSION BP(0:JCAP,0:JCAP)
      DIMENSION AQR(0:JCAP,0:JCAP)
      DIMENSION BQR(0:JCAP,0:JCAP)
      DIMENSION GR(0:JCAP,0:JCAP)
      DIMENSION DEL2(0:JCAP,0:JCAP)
      INTEGER J1X, J2X, J3X
C--------
      RERTH=CONMC('RERTH$')
      DO J2X = 1, JCAP*(JCAP + 2) + 1
         AP(J2X-1,0) = 0.
         BP(J2X-1,0) = 0.
         AQR(J2X-1,0) = 0.
         BQR(J2X-1,0) = 0.
         GR(J2X-1,0) = 0.
         DEL2(J2X-1,0) = 0.
      END DO
      DO 20 M=0,JCAP
        DO 10 L=0,JCAP-M
          N=M+L
          AP(M,L)=SQRT((2.*N+1.)*(2.*N+3.)/
     *         ((N-L+1.)*(N+L+1.)))
          BP(M,L)=-SQRT((N-L)*(N+L)*(2.*N+3.)/
     *         ((N-L+1.)*(N+L+1.)*MAX(1.,2.*N-1.)))
          AQR(M,L)=AP(M,L)*N/(N+2.)
          BQR(M,L)=BP(M,L)*N*(N-1.)/((N+1.)*(N+2.))
          GR(M,L)=RERTH*AP(M,L)/((N+1.)*(N+2.))
          DEL2(M,L)=N*(N+1.)/(RERTH*RERTH)
10      CONTINUE
20    CONTINUE
      DO 40 M=0,JCAP-1
        MII=JCAP-M
      IF (M.NE.MII .OR. JCAP-M.LT.1+M .OR. JCAP.LT.1) THEN
         DO 30 L = 1, JCAP - M
            DEL2(MII,JCAP+1-L) = DEL2(M,L)
   30    CONTINUE
      ELSE
         DO L = 1, JCAP - M
            LII = JCAP + 1 - L
            DEL2(MII,LII) = DEL2(M,L)
         END DO
      ENDIF
40    CONTINUE
      RETURN
      END
