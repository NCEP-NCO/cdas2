      SUBROUTINE HTOPER(ZS,DS,HS,QS,PS,U,V,VORT,T,P,PLON,PLAT,Q,
     *     BHALF,BHALFP,NSIG,JCAP,NLON,NLATH,DEL2,
     *     PLN,QLN,RLN,TRIGS,IFAX,
     *     AGVZ,WGVZ,BVZ,NMDSZH,VZ,VD,VH,VQ,IN,BALN)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    HTOPER      TRANSPOSE OF HOPER
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-06
C
C ABSTRACT: APPLY TRANSPOSE OF HOPER, GOING FROM GRID TO SPECTRAL.
C
C PROGRAM HISTORY LOG:
C   90-10-06  PARRISH
C
C   INPUT ARGUMENT LIST:
C     U,V,VORT,T,P,PLON,PLAT,Q - U,V,ETC ON GAUSSIAN GRID
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
C
C   OUTPUT ARGUMENT LIST:
C     ZS,DS,HS,QS,PS - COEFS OF VORT, DIV, UNBAL T, UNBAL LOG(PS), Q
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION AGVZ(0:JCAP,NSIG,NMDSZH)
C-CRA          DIMENSION WGVZ(0:JCAP,NMDSZH)
C-CRA          DIMENSION BVZ(0:JCAP,NSIG,NMDSZH)
C-CRA          DIMENSION BHALF((JCAP+1)*(JCAP+2),NSIG,4)
C-CRA          DIMENSION BHALFP((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION ZS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION DS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION HS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION QS((JCAP+1)*(JCAP+2),NSIG)
C-CRA          DIMENSION PS((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION U(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION V(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION VORT(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION Q(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION VZ(NSIG,NSIG),VD(NSIG,NSIG)
C-CRA          DIMENSION VQ(NSIG,NSIG),VH(NSIG,NSIG)
C-CRA          DIMENSION DEL2((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION TRIGS(NLON*2),IFAX(10)
C-CRA          DIMENSION PLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA          DIMENSION IN((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION BALN((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION PSD((JCAP+1)*(JCAP+2))
C-CRA          DIMENSION ZSF((JCAP+1)*(JCAP+2),NMDSZH)
C-CRA          DIMENSION WORK((JCAP+1)*(JCAP+2),NSIG)
C-CRA          REAL T(2*NLATH+1,NLON+2,NSIG),P(2*NLATH+1,NLON+2)
C-CRA          REAL PLON(2*NLATH+1,NLON+2),PLAT(2*NLATH+1,NLON+2)
 
          DIMENSION AGVZ(0:62,28,28)
          DIMENSION WGVZ(0:62,28)
          DIMENSION BVZ(0:62,28,28)
          DIMENSION BHALF((62+1)*(62+2),28,4)
          DIMENSION BHALFP((62+1)*(62+2))
          DIMENSION ZS((62+1)*(62+2),28)
          DIMENSION DS((62+1)*(62+2),28)
          DIMENSION HS((62+1)*(62+2),28)
          DIMENSION QS((62+1)*(62+2),28)
          DIMENSION PS((62+1)*(62+2))
          DIMENSION U(2*48+1,192+2,28)
          DIMENSION V(2*48+1,192+2,28)
          DIMENSION VORT(2*48+1,192+2,28)
          DIMENSION Q(2*48+1,192+2,28)
          DIMENSION VZ(28,28),VD(28,28)
          DIMENSION VQ(28,28),VH(28,28)
          DIMENSION DEL2((62+1)*(62+2))
          DIMENSION TRIGS(192*2),IFAX(10)
          DIMENSION PLN((62+1)*(62+2),48)
          DIMENSION QLN((62+1)*(62+2),48)
          DIMENSION RLN((62+1)*(62+2),48)
          DIMENSION IN((62+1)*(62+2))
          DIMENSION BALN((62+1)*(62+2))
          DIMENSION PSD((62+1)*(62+2))
          DIMENSION ZSF((62+1)*(62+2),28)
          DIMENSION WORK((62+1)*(62+2),28)
          REAL T(2*48+1,192+2,28),P(2*48+1,192+2)
          REAL PLON(2*48+1,192+2),PLAT(2*48+1,192+2)
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C--------
         NC=(JCAP+1)*(JCAP+2)
         NG=(2*NLATH+1)*(NLON+2)
         DO KK=1,NSIG*3+2
          IF(KK.EQ.3*NSIG+1)
     *      CALL TS2GRAD(PSD,PLON,PLAT,JCAP,NLON,NLATH,QLN,RLN,
     *            TRIGS,IFAX)
          IF(KK.EQ.3*NSIG+2)
     *      CALL TS2G0(PS,P,JCAP,NLON,NLATH,PLN,TRIGS,IFAX)
          K=MOD(KK-1,NSIG)+1
          IF(KK.GE.1.AND.KK.LE.NSIG) THEN
           CALL TS2G0(ZS(1,K),VORT(1,1,K),JCAP,NLON,NLATH,PLN,
     *             TRIGS,IFAX)
           CALL TS2GVEC(WORK(1,K),DS(1,K),U(1,1,K),V(1,1,K),
     *          JCAP,NLON,NLATH,QLN,RLN,TRIGS,IFAX)
           DO I=1,NC
            ZS(I,K)=ZS(I,K)+WORK(I,K)
           END DO
          END IF
          IF(KK.GE.NSIG+1.AND.KK.LE.2*NSIG)
     *      CALL TS2G0(HS(1,K),T(1,1,K),JCAP,NLON,NLATH,PLN,
     *              TRIGS,IFAX)
          IF(KK.GE.2*NSIG+1.AND.KK.LE.3*NSIG)
     *      CALL TS2G0(QS(1,K),Q(1,1,K),JCAP,NLON,NLATH,PLN,
     *              TRIGS,IFAX)
         END DO
C-CRA                   PS=PS-DEL2*PSD
C       DIMENSION PS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          PS(ITMP)=PS(ITMP)-DEL2(ITMP)*PSD(ITMP)
          ENDDO
C--------
C-------- NEXT DO VERTICAL TRANSFORMS
C--------
C-------- TSUM IN VERTICAL ZS                       
      DO J=1,NSIG
         DO I=1,NC
          WORK(I,J)=VZ(1,J)*ZS(I,1)
         END DO
        DO K=2,NSIG
         DO I=1,NC
          WORK(I,J)=VZ(K,J)*ZS(I,K)
     *             +WORK(I,J)
         END DO
        END DO
      END DO
C--------
C-------- TSUM IN VERTICAL QS                       
      DO J=1,NSIG
         DO I=1,NC
          ZS(I,J)=WORK(I,J)
         END DO
         DO I=1,NC
          WORK(I,J)=VQ(1,J)*QS(I,1)
         END DO
        DO K=2,NSIG
         DO I=1,NC
          WORK(I,J)=VQ(K,J)*QS(I,K)
     *             +WORK(I,J)
         END DO
        END DO
      END DO
C-CRA                QS=WORK
C       DIMENSION QS((JCAP+1)*(JCAP+2),NSIG)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG
          QS(ITMP,1)=WORK(ITMP,1)
          ENDDO
C-CRA                WORK=DS
C       DIMENSION WORK((JCAP+1)*(JCAP+2),NSIG)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG
          WORK(ITMP,1)=DS(ITMP,1)
          ENDDO
C-CRA                DS=0.
C       DIMENSION DS((JCAP+1)*(JCAP+2),NSIG)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG
          DS(ITMP,1)=0.
          ENDDO
C-CRA                ZSF=0.
C       DIMENSION ZSF((JCAP+1)*(JCAP+2),NMDSZH)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NMDSZH
          ZSF(ITMP,1)=0.
          ENDDO
C--------------- TSUM IN VERTICAL DS                       
         DO J=1,NSIG
          DO K=1,NSIG
           DO I=1,NC
            DS(I,J)=DS(I,J)+WORK(I,K)*VD(K,J)
           END DO
          END DO
          IF(J.LE.NMDSZH) THEN
           DO K=1,NSIG
            DO I=1,NC
             ZSF(I,J)=ZSF(I,J)+BVZ(IN(I),K,J)*WORK(I,K)
            END DO
           END DO
          END IF
         END DO
C---------------DO TTEMP AND TPSFC
C-CRA                   WORK=HS
C       DIMENSION WORK((JCAP+1)*(JCAP+2),NSIG)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG
          WORK(ITMP,1)=HS(ITMP,1)
          ENDDO
C-CRA                   HS=0.
C       DIMENSION HS((JCAP+1)*(JCAP+2),NSIG)
          DO ITMP=1,((JCAP+1)*(JCAP+2))*NSIG
          HS(ITMP,1)=0.
          ENDDO
         DO J=1,NSIG
          IF(J.LE.NMDSZH) THEN
           DO I=1,NC
            ZSF(I,J)=ZSF(I,J)+WGVZ(IN(I),J)*PS(I)
           END DO
           DO K=1,NSIG
            DO I=1,NC
             ZSF(I,J)=ZSF(I,J)+AGVZ(IN(I),K,J)*WORK(I,K)
            END DO
           END DO
          END IF
          DO K=1,NSIG
           DO I=1,NC
            HS(I,J)=HS(I,J)+WORK(I,K)*VH(K,J)
           END DO
          END DO
         END DO
C------------------------TAPPLY SPECTRAL BALANCE OPERATOR
C------------------------TO ZS
        DO K=1,NMDSZH
         II0=2*(JCAP+1)
         IM0=0
         DO M=1,JCAP
          DO LL=1,2*(JCAP+1-M)
           ZS(IM0+LL,K)=ZS(IM0+LL,K)+BALN(II0+LL)*ZSF(II0+LL,K)
          END DO
          II0=II0+2*(JCAP+1-M)
          IM0=IM0+2*(JCAP+2-M)
         END DO
         II0=0
         IP0=2*(JCAP+1)
         DO M=0,JCAP-1
          DO LL=1,2*(JCAP-M)
           ZS(IP0+LL,K)=ZS(IP0+LL,K)+BALN(IP0+LL)*ZSF(II0+LL,K)
          END DO
          II0=II0+2*(JCAP+1-M)
          IP0=IP0+2*(JCAP-M)
         END DO
        END DO
      DO J=1,NSIG
       IF(J .EQ. 1)THEN
         DO I=1,NC
           PS(I)=PS(I)*BHALFP(I)
         END DO
       END IF
       DO I=1,NC
        ZS(I,J)=ZS(I,J)*BHALF(I,J,1)
        DS(I,J)=DS(I,J)*BHALF(I,J,2)
        HS(I,J)=HS(I,J)*BHALF(I,J,3)
        QS(I,J)=QS(I,J)*BHALF(I,J,4)
       END DO
      END DO
      RETURN
      END
