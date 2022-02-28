       SUBROUTINE GETBALN(BALN,JCAP)
C--------
C-------- OBTAIN CONSTANTS FOR SPECTRAL APPLICATION OF LINEAR
C-------- BALANCE OPERATOR--- DEL**(-2) DEL DOT F DEL DEL**(-2)
C--------
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    GETBALN    GET CONSTS FOR SPECTRAL LIN-BAL OPERATOR.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 94-02-11
C
C ABSTRACT: GET CONSTS FOR SPECTRAL LINEAR-BALANCE OPERATOR.
C
C PROGRAM HISTORY LOG:
C   94-02-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     JCAP     - TRIANGULAR TRUNCATION
C
C   OUTPUT ARGUMENT LIST:
C     BALN     - BALANCE OPERATOR CONSTANTS.
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA             DIMENSION BALN((JCAP+1)*(JCAP+2))
 
             DIMENSION BALN((62+1)*(62+2))
C--------
         RERTH=CONMC('RERTH$')
         ABC=-2.*CONMC('OMEGA$')*RERTH**2/CONMC('G$')
C-CRA                   BALN(1:2*(JCAP+1))=0.
          DO ITMP=1,2*(JCAP+1)
                   BALN(ITMP)=0.
          END DO
         II=2*(JCAP+1)-1
         DO M=1,JCAP
          II=II+2
          RN=M
          EPS=SQRT(RN**2/(4.*RN**2-1.))
          BALN(II)=ABC*EPS/(RN*RN)
          BALN(II+1)=0.
          IF(M.LT.JCAP) THEN
           DO L=1,JCAP-M
            II=II+2
            RL=L
            RN=M+L
            EPS=SQRT((RN**2-RL**2)/(4.*RN**2-1.))
            BALN(II)=ABC*EPS/(RN*RN)
            BALN(II+1)=BALN(II)
           END DO
          END IF
         END DO
         BALN(2*(JCAP+1)+1)=0.
       RETURN
       END
