       SUBROUTINE HOPER(ZS,DS,HS,QS,PS,U,V,VORT,T,P,PLON,PLAT,Q,
     *     BHALF,BHALFP,NSIG,JCAP,NLON,NLATH,DEL2,
     *     PLN,QLN,RLN,TRIGS,IFAX,
     *     AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    HOPER      ANALYSIS VARIABLES TO GRID VARIABLES
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-06
C
C ABSTRACT: CONVERT ANALYSIS VARIABLES TO GRID VARIABLES
C
C PROGRAM HISTORY LOG:
C   90-10-06  PARRISH
C   94-02-02  PARRISH
C
C   INPUT ARGUMENT LIST:
C     ZS,DS,HS,QS,PS - COEFS OF VORT, DIV, UNBAL T, UNBAL LOG(PS), Q
C     BHALF    - BACKGROUND ERROR STATS
C     BHALFP   - BACKGROUND ERROR STATS SURFACE PRESSURE
C     NSIG     - NUMBER OF SIGMA LEVELS
C     JCAP     - TRIANGULAR TRUNCATION
C     NLON     - NUMBER OF LONGITUDES
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     DEL2     - N*(N+1)/A**2
C     TRIGS,IFAX - USED BY FFT
C     AGVZ     - MASS-VARIABLE MODES TO TEMPERATURE CONVERSION
C     WGVZ     - MASS-VARIABLE MODES TO LOG(PSFC) CONVERSION
C     BVZ      - MASS-VARIABLE MODES TO DIVERGENCE CONVERSION
C     NMDSZH   - NUMBER OF MODES USED IN BALANCE EQN.
C     VZ       - VERTICAL MODE MATRIX - Z    
C     VD       - VERTICAL MODE MATRIX - D    
C     VH       - VERTICAL MODE MATRIX - TEMPS
C     VQ       - VERTICAL MODE MATRIX - Q   
C     IN       - TOTAL WAVENUMBER INDEX ARRAY
C     BALN     - SPECTRAL BALANCE OPERATOR CONSTANTS
C
C   OUTPUT ARGUMENT LIST:
C     U,V,VORT,T,P,PLON,PLAT,Q - U,V,ETC ON GAUSSIAN GRID
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C-CRA             REAL T(2*NLATH+1,NLON+2,NSIG),P(2*NLATH+1,NLON+2)
C-CRA             REAL PLON(2*NLATH+1,NLON+2),PLAT(2*NLATH+1,NLON+2)
C-CRA             DIMENSION AGVZ(0:JCAP,NSIG,NMDSZH)
C-CRA             DIMENSION WGVZ(0:JCAP,NMDSZH)
C-CRA             DIMENSION BVZ(0:JCAP,NSIG,NMDSZH)
C-CRA             DIMENSION U(2*NLATH+1,NLON+2,NSIG)
C-CRA             DIMENSION V(2*NLATH+1,NLON+2,NSIG)
C-CRA             DIMENSION VORT(2*NLATH+1,NLON+2,NSIG)
C-CRA             DIMENSION Q(2*NLATH+1,NLON+2,NSIG)
C-CRA             DIMENSION VZ(NSIG,NSIG),VD(NSIG,NSIG)
C-CRA             DIMENSION VQ(NSIG,NSIG),VH(NSIG,NSIG)
C-CRA             DIMENSION BHALF((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA             DIMENSION BHALFP((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION ZS((JCAP+1)*(JCAP+2),NSIG)
C-CRA             DIMENSION DS((JCAP+1)*(JCAP+2),NSIG)
C-CRA             DIMENSION HS((JCAP+1)*(JCAP+2),NSIG)
C-CRA             DIMENSION QS((JCAP+1)*(JCAP+2),NSIG)
C-CRA             DIMENSION PS((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION DEL2((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION TRIGS(NLON*2),IFAX(10)
C-CRA             DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION IN((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION BALN((JCAP+1)*(JCAP+2))
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C-CRA             DIMENSION PSD((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION ZSF((JCAP+1)*(JCAP+2),NMDSZH)
C-CRA             DIMENSION WORK((JCAP+1)*(JCAP+2),NSIG)

             REAL T(2*48+1,192+2,28),P(2*48+1,192+2)
             REAL PLON(2*48+1,192+2),PLAT(2*48+1,192+2)
             DIMENSION AGVZ(0:62,28,28)
             DIMENSION WGVZ(0:62,28)
             DIMENSION BVZ(0:62,28,28)
             DIMENSION U(2*48+1,192+2,28)
             DIMENSION V(2*48+1,192+2,28)
             DIMENSION VORT(2*48+1,192+2,28)
             DIMENSION Q(2*48+1,192+2,28)
             DIMENSION VZ(28,28),VD(28,28)
             DIMENSION VQ(28,28),VH(28,28)
             DIMENSION BHALF((62+1)*(62+2),28,4)
             DIMENSION BHALFP((62+1)*(62+2))
             DIMENSION ZS((62+1)*(62+2),28)
             DIMENSION DS((62+1)*(62+2),28)
             DIMENSION HS((62+1)*(62+2),28)
             DIMENSION QS((62+1)*(62+2),28)
             DIMENSION PS((62+1)*(62+2))
             DIMENSION DEL2((62+1)*(62+2))
             DIMENSION TRIGS(192*2),IFAX(10)
             DIMENSION PLN((62+1)*(62+2),48)
             DIMENSION QLN((62+1)*(62+2),48)
             DIMENSION RLN((62+1)*(62+2),48)
             DIMENSION IN((62+1)*(62+2))
             DIMENSION BALN((62+1)*(62+2))
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
             DIMENSION PSD((62+1)*(62+2))
             DIMENSION ZSF((62+1)*(62+2),28)
             DIMENSION WORK((62+1)*(62+2),28)
C--------
         NC=(JCAP+1)*(JCAP+2)
         NG=(2*NLATH+1)*(NLON+2)
C--------
C-------- FIRST SUM IN VERTICAL, AND ZERO VARIOUS ARRAYS)
         DO K=1,NSIG
          IF(K .EQ. 1)THEN
C-CRA               P=0.
C-CRA               PLON=0.
C-CRA               PLAT=0.
               DO I=1,(2*NLATH+1)*(NLON+2)
               P(I,1)=0.
               PLON(I,1)=0.
               PLAT(I,1)=0.
               ENDDO
           DO I=1,NC
            PS(I)=PS(I)*BHALFP(I)
           END DO
          END IF
          DO I=1,NC
           ZS(I,K)=ZS(I,K)*BHALF(I,K,1)
           DS(I,K)=DS(I,K)*BHALF(I,K,2)
           HS(I,K)=HS(I,K)*BHALF(I,K,3)
           QS(I,K)=QS(I,K)*BHALF(I,K,4)
          END DO
          DO I=1,NG
           U(I,1,K)=0.
           V(I,1,K)=0.
           VORT(I,1,K)=0.
           T(I,1,K)=0.
           Q(I,1,K)=0.
          END DO
         END DO
C------------------------APPLY SPECTRAL BALANCE OPERATOR
C------------------------TO ZS
C-CRA            ZSF=0.
            DO I=1,(JCAP+1)*(JCAP+2)*NMDSZH
            ZSF(I,1)=0.
            ENDDO
        DO K=1,NMDSZH
         II0=2*(JCAP+1)
         IM0=0
         DO M=1,JCAP
          DO LL=1,2*(JCAP+1-M)
           ZSF(II0+LL,K)=ZSF(II0+LL,K)+BALN(II0+LL)*ZS(IM0+LL,K)
          END DO
          II0=II0+2*(JCAP+1-M)
          IM0=IM0+2*(JCAP+2-M)
         END DO
         II0=0
         IP0=2*(JCAP+1)
         DO M=0,JCAP-1
          DO LL=1,2*(JCAP-M)
           ZSF(II0+LL,K)=ZSF(II0+LL,K)+BALN(IP0+LL)*ZS(IP0+LL,K)
          END DO
          II0=II0+2*(JCAP+1-M)
          IP0=IP0+2*(JCAP-M)
         END DO
        END DO
C---------------DO TEMP AND PSFC
C-CRA             WORK=0.
             DO I=1,(JCAP+1)*(JCAP+2)*NSIG
             WORK(I,1)=0.
             ENDDO
         DO K=1,NSIG
          IF(K .EQ. 1)THEN
           DO J=1,NMDSZH
            DO I=1,NC
             PS(I)=PS(I)
     *             +WGVZ(IN(I),J)*ZSF(I,J)
            END DO
           END DO
          END IF
          DO J=1,NMDSZH
           DO I=1,NC
            WORK(I,K)=WORK(I,K)
     *             +AGVZ(IN(I),K,J)*ZSF(I,J)
           END DO
          END DO
          DO J=1,NSIG
           DO I=1,NC
            WORK(I,K)=WORK(I,K)+VH(K,J)*HS(I,J)
           END DO
          END DO
         END DO
         DO I=1,NSIG*NC
          HS(I,1)=WORK(I,1)
          WORK(I,1)=0.
         END DO
C--------------- SUM IN VERTICAL DS                       
         DO K=1,NSIG
          DO J=1,NSIG
           DO I=1,NC
            WORK(I,K)=WORK(I,K)+VD(K,J)*DS(I,J)
           END DO
          END DO
          DO J=1,NMDSZH
           DO I=1,NC
            WORK(I,K)=WORK(I,K)+BVZ(IN(I),K,J)*ZSF(I,J)
           END DO
          END DO
         END DO
         DO I=1,NC*NSIG
          DS(I,1)=WORK(I,1)
          WORK(I,1)=0.
         END DO
C--------
C-------- SUM IN VERTICAL QS                       
         DO K=1,NSIG
          DO J=1,NSIG
           DO I=1,NC
            WORK(I,K)=WORK(I,K)+VQ(K,J)*QS(I,J)
           END DO
          END DO
         END DO
         DO I=1,NSIG*NC
          QS(I,1)=WORK(I,1)
          WORK(I,1)=0.
         END DO
C-------
C-------- SUM IN VERTICAL ZS                       
         DO K=1,NSIG
          DO J=1,NSIG
           DO I=1,NC
            WORK(I,K)=WORK(I,K)+VZ(K,J)*ZS(I,J)
           END DO
          END DO
         END DO
         DO I=1,NSIG*NC
          ZS(I,1)=WORK(I,1)
         END DO
         DO I=1,NC
          PSD(I)=-DEL2(I)*PS(I)
         END DO
         DO KK=1,NSIG*3+2
          IF(KK.EQ.3*NSIG+1)
     *      CALL S2GRAD(PSD,PLON,PLAT,JCAP,NLON,NLATH,QLN,RLN,
     *           TRIGS,IFAX)
          IF(KK.EQ.3*NSIG+2)
     *      CALL S2G0(PS,P,JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
          K=MOD(KK-1,NSIG)+1
          IF(KK.GE.1.AND.KK.LE.NSIG) THEN
           CALL S2G0(ZS(1,K),VORT(1,1,K),JCAP,NLON,NLATH,PLN,
     *              TRIGS,IFAX)
           CALL S2GVEC(ZS(1,K),DS(1,K),U(1,1,K),V(1,1,K),
     *        JCAP,NLON,NLATH,QLN,RLN,TRIGS,IFAX)
          END IF
          IF(KK.GE.NSIG+1.AND.KK.LE.2*NSIG)
     *      CALL S2G0(HS(1,K),T(1,1,K),JCAP,NLON,NLATH,PLN,
     *              TRIGS,IFAX)
          IF(KK.GE.2*NSIG+1.AND.KK.LE.3*NSIG)
     *      CALL S2G0(QS(1,K),Q(1,1,K),JCAP,NLON,NLATH,PLN,
     *              TRIGS,IFAX)
         END DO
       RETURN
       END
