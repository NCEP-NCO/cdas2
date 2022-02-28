       SUBROUTINE GTBDIVT(IDIVT,BDIVT,NSIG,JCAP,NSIGDIVT,
     *            JCAPDIVT,SIGL)
C-------------
C-------------BRING IN ESTIMATES OF DIVTEND ERROR VARIANCE,
C-----------  INTERPOLATE IN VERTICAL, AND TRUNCATE IN HORIZONTAL
C-----------
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    GTBDIVT    READ DIVTEND ERROR VARIANCES.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 94-02-11
C
C ABSTRACT: READ DIVTEND ERROR VARIANCE, INTERPOLATE IN VERT, TRUNCATE
C            IN HORIZONTAL, AS NECESSARY.
C
C PROGRAM HISTORY LOG:
C   94-02-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     IDIVT    - INPUT UNIT NUMBER CONTAINING DIVTEND ERROR VARIANCES
C     NSIG     - NUMBER OF SIGMA LEVELS
C     JCAP     - TRIANGULAR TRUNCATION
C     NSIGDIVT - NUMBER OF SIGMA LEVELS IN INPUT DIVTEND ERRORS
C     JCAPDIVT - TRIANGULAR TRUNCATION FOR INPUT DIVTEND ERRORS
C     SIGL     - ANALYSIS SIGMA LEVELS
C
C   OUTPUT ARGUMENT LIST:
C     BDIVT    - INVERSE OF DIVTEND ERROR VARIANCES.
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C-CRA             DIMENSION BDIVT(0:JCAP,NSIG),SIGL(NSIG)
C-CRA             DIMENSION CDIVT(0:JCAPDIVT,NSIGDIVT),SIGLDIVT(NSIGDIVT)
C-CRA             DIMENSION TDIVT(0:JCAPDIVT,NSIG)
C-CRA             DIMENSION RLSG(NSIG),RLSGDIVT(NSIGDIVT)
C-CRA             DIMENSION GRID(NSIG)

             DIMENSION BDIVT(0:62,28),SIGL(28)
             DIMENSION CDIVT(0:126,28)
             DIMENSION SIGLDIVT(28)
             DIMENSION TDIVT(0:126,28)
             DIMENSION RLSG(28),RLSGDIVT(28)
             DIMENSION GRID(28)
C-----------
         REWIND IDIVT
         READ(IDIVT)MSIGDIVT,MCAPDIVT,SIGLDIVT,CDIVT
         CLOSE(IDIVT)
C-CRA             RLSGDIVT=LOG(SIGLDIVT)
             DO I=1,NSIGDIVT
             RLSGDIVT(I)=LOG(SIGLDIVT(I))
             ENDDO
C-CRA             RLSG=LOG(SIGL)
             DO I=1,NSIG
             RLSG(I)=LOG(SIGL(I))
             ENDDO
C-CRA             GRID=RLSG
             DO I=1,NSIG
             GRID(I)=RLSG(I)
             ENDDO
         CALL GDCRDN(GRID,NSIG,RLSGDIVT,NSIGDIVT)
         DO K=1,NSIG
          I0=GRID(K)
          I0=MAX(1,MIN(NSIGDIVT-1,I0))
          I1=I0+1
          DEL=GRID(K)-I0
          W0=1.-DEL
          W1=DEL
          DO N=0,JCAPDIVT
           TDIVT(N,K)=W0*CDIVT(N,I0)+W1*CDIVT(N,I1)
          END DO
          DO N=0,MIN(JCAP,JCAPDIVT)
           BDIVT(N,K)=TDIVT(N,K)
          END DO
          IF(JCAP.GT.JCAPDIVT) THEN
           DO N=JCAPDIVT+1,JCAP
            BDIVT(N,K)=TDIVT(JCAPDIVT,K)
           END DO
          END IF
          DO N=1,JCAP
           BDIVT(N,K)=1./BDIVT(N,K)
          END DO
         END DO
       RETURN
       END