      SUBROUTINE RDGESC(ZC,DC,TC,QC,PC,RC,HOURG,IDATEG,SIGI,SIGL,
     *  INGES,JCAP,NSIG,ON85,ML2LM,FACTSLM,FACTVLM)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    RDGESC     READ SIGMA COEFS AND REORDER.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-10
C
C ABSTRACT: READ GUESS SIGMA COEFS, AND REORDER TO INTERNAL FORMAT.
C
C PROGRAM HISTORY LOG:
C   90-10-10  PARRISH
C
C   INPUT ARGUMENT LIST:
C     INGES    - UNIT NUMBER OF GUESS COEFS
C     JCAP     - TRIANGULAR TRUNCATION
C     NSIG     - NUMBER OF SIGMA LEVELS
C
C   OUTPUT ARGUMENT LIST:
C     ZC,DC,TC,QC,PC,RC - GES SIG COEFS OF VORT,DIV,T,Q,LN(PS),Z0
C     HOURG    - GUESS FORECAST HOUR
C     IDATEG   - INITIAL DATE OF GUESS
C     SIGI     - SIGMA VALUES AT INTERFACE OF EACH SIGMA LAYER
C     SIGL     - SIGMA VALUES AT MID-POINT OF EACH SIGMA LAYER
C     ON85     - ON85 DATE RECORD FOR GUESS COEFS
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION ZC((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION DC((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION TC((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION QC((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION PC((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION RC((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION IDATEG(4),SIGI(NSIG+1),SIGL(NSIG)
C-CRA          DIMENSION ML2LM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTSLM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION FACTVLM((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION Z((JCAP+1)*(JCAP+2))
C-CRA          CHARACTER*4 ON85(8)
 
          DIMENSION ZC((62+1)*(62+2),28)
          DIMENSION DC((62+1)*(62+2),28)
          DIMENSION TC((62+1)*(62+2),28)
          DIMENSION QC((62+1)*(62+2),28)
          DIMENSION PC((62+1)*(62+2))
          DIMENSION RC((62+1)*(62+2))
          DIMENSION IDATEG(4),SIGI(28+1),SIGL(28)
          DIMENSION ML2LM((62+1)*(62+2))
          DIMENSION FACTSLM((62+1)*(62+2))
          DIMENSION FACTVLM((62+1)*(62+2))
          DIMENSION Z((62+1)*(62+2))
          CHARACTER*4 ON85(8)
C--------
C-------- LOCAL SPACE
C--------
C-------
      NC=(JCAP+1)*(JCAP+2)
      REWIND INGES
C-------- HOUR,IDATE, ETC.
      READ(INGES)ON85
      READ(INGES)HOURG,IDATEG,SIGI,SIGL
C-------- TERRAIN COEFS
      READ(INGES)Z
      DO I=1,NC
       RC(I)=FACTSLM(I)*Z(ML2LM(I))
      END DO
C-------- SFCP COEFFICIENTS
      READ(INGES)Z
      DO I=1,NC
       PC(I)=FACTSLM(I)*Z(ML2LM(I))
      END DO
C-------- TEMP COEFFICIENTS
      DO K=1,NSIG
        READ(INGES)Z
        DO I=1,NC
         TC(I,K)=FACTSLM(I)*Z(ML2LM(I))
        END DO
      END DO
C-------- DIV AND VORT
      DO K=1,NSIG
        READ(INGES)Z
        DO I=1,NC
         DC(I,K)=FACTVLM(I)*Z(ML2LM(I))
        END DO
        READ(INGES)Z
        DO I=1,NC
         ZC(I,K)=FACTVLM(I)*Z(ML2LM(I))
        END DO
      END DO
C-------- Q COEFS
      DO K=1,NSIG
        READ(INGES)Z
        DO I=1,NC
         QC(I,K)=FACTSLM(I)*Z(ML2LM(I))
        END DO
      END DO
      WRITE(6,700)JCAP,NSIG,HOURG,IDATEG
700   FORMAT(' GUESS SIGMA COEFFICIENTS READ IN, JCAP,NSIG=',
     *  2I6,/,' HOUR,IDATE=',F10.1,4I4)
      CLOSE (INGES)
      RETURN
      END