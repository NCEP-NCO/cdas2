       SUBROUTINE S2GVEC(ZS,DS,U,V,JCAP,NLON,NLATH,QLN,RLN,TRIGS,IFAX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    S2GVEC     VORT, DIV TO U,V
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 94-04-08
C
C ABSTRACT: SPECTRAL VORT, DIV COEFS TO GRID U,V.
C
C PROGRAM HISTORY LOG:
C   94-04-08  PARRISH
C
C   INPUT ARGUMENT LIST:
C     ZS       - VORTICITY COEFFICIENTS
C     DS       - DIVERGENCE COEFFICIENTS
C     JCAP     - TRIANGULAR TRUNCATION
C     NLON     - NUMBER OF LONGITUDES
C     NLATH    - NUMBER OF GAUSSIAN LATS IN ONE HEMISPHERE
C     QLN      - Q(N,L)
C     RLN      - R(N,L)
C
C   OUTPUT ARGUMENT LIST:
C     U        - LONGITUDE COMPONENT OF WINDS
C     V        - LATITUDE COMPONENT OF WINDS
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA             DIMENSION ZS((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION DS((JCAP+1)*(JCAP+2))
C-CRA             DIMENSION U(2*NLATH+1,NLON+2)
C-CRA             DIMENSION V(2*NLATH+1,NLON+2)
C-CRA             DIMENSION QLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION RLN((JCAP+1)*(JCAP+2),NLATH)
C-CRA             DIMENSION TRIGS(2*NLON),IFAX(10)
C-CRA             DIMENSION WORK(2*(2*NLATH+1)*(NLON+2))
C-CRA             DIMENSION UE(2*(JCAP+1)),UO(2*(JCAP+1))
C-CRA             DIMENSION VE(2*(JCAP+1)),VO(2*(JCAP+1))
 
             DIMENSION ZS((62+1)*(62+2))
             DIMENSION DS((62+1)*(62+2))
             DIMENSION U(2*48+1,192+2)
             DIMENSION V(2*48+1,192+2)
             DIMENSION QLN((62+1)*(62+2),48)
             DIMENSION RLN((62+1)*(62+2),48)
             DIMENSION TRIGS(2*192),IFAX(10)
             DIMENSION WORK(2*(2*48+1)*(192+2))
             DIMENSION UE(2*(62+1)),UO(2*(62+1))
             DIMENSION VE(2*(62+1)),VO(2*(62+1))
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C---------------
C-CRA                   U=0.
C          DIMENSION U(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          U(ITMP,1)=0.
          ENDDO
C-CRA                   V=0.
C          DIMENSION V(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          V(ITMP,1)=0.
          ENDDO
         DO J=1,NLATH
          JR=2*NLATH+1-J
          II0=0
C-CRA                    UE=0.
C          DIMENSION UE(2*(JCAP+1)),UO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          UE(ITMP)=0.
          ENDDO
C-CRA                    UO=0.
C          DIMENSION UE(2*(JCAP+1)),UO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          UO(ITMP)=0.
          ENDDO
C-CRA                    VE=0.
C          DIMENSION VE(2*(JCAP+1)),VO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          VE(ITMP)=0.
          ENDDO
C-CRA                    VO=0.
C          DIMENSION VE(2*(JCAP+1)),VO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          VO(ITMP)=0.
          ENDDO
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            UO(LL)=UO(LL)+RLN(II0+LL,J)*ZS(II0+LL)
            VO(LL)=VO(LL)-RLN(II0+LL,J)*DS(II0+LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             UE(LL)=UE(LL)+RLN(II0+LL,J)*ZS(II0+LL)
             VE(LL)=VE(LL)-RLN(II0+LL,J)*DS(II0+LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
C----------
C---------- NOW COMBINE EVEN AND ODD PARTS
C----------
          DO LL=1,2*(JCAP+1)
           U(J,LL)=UE(LL)+UO(LL)
           U(JR,LL)=UE(LL)-UO(LL)
           V(J,LL)=VE(LL)+VO(LL)
           V(JR,LL)=VE(LL)-VO(LL)
          END DO
         END DO
C---------------
C-------------MULTIPLY DIV,VORT BY I
         DO I=1,(JCAP+1)*(JCAP+2),2
          DIVR=DS(I)
          DIVI=DS(I+1)
          DS(I)=DIVI
          DS(I+1)=-DIVR
          VORR=ZS(I)
          VORI=ZS(I+1)
          ZS(I)=VORI
          ZS(I+1)=-VORR
         END DO
         DO J=1,NLATH
          JR=2*NLATH+1-J
          II0=0
C-CRA                    UE=0.
C          DIMENSION UE(2*(JCAP+1)),UO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          UE(ITMP)=0.
          ENDDO
C-CRA                    UO=0.
C          DIMENSION UE(2*(JCAP+1)),UO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          UO(ITMP)=0.
          ENDDO
C-CRA                    VE=0.
C          DIMENSION VE(2*(JCAP+1)),VO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          VE(ITMP)=0.
          ENDDO
C-CRA                    VO=0.
C          DIMENSION VE(2*(JCAP+1)),VO(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          VO(ITMP)=0.
          ENDDO
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            UE(LL)=UE(LL)+QLN(II0+LL,J)*DS(II0+LL)
            VE(LL)=VE(LL)+QLN(II0+LL,J)*ZS(II0+LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             UO(LL)=UO(LL)+QLN(II0+LL,J)*DS(II0+LL)
             VO(LL)=VO(LL)+QLN(II0+LL,J)*ZS(II0+LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
C----------
C---------- NOW COMBINE EVEN AND ODD PARTS
C----------
          DO LL=1,2*(JCAP+1)
           U(J,LL)=U(J,LL)+UE(LL)+UO(LL)
           U(JR,LL)=U(JR,LL)+UE(LL)-UO(LL)
           V(J,LL)=V(J,LL)+VE(LL)+VO(LL)
           V(JR,LL)=V(JR,LL)+VE(LL)-VO(LL)
          END DO
         END DO
C--------
C-------- FINALLY DO FOURIER SUMS IN LONGITUDE
C--------
         LOT=NLATH*2
         NLAX=LOT+1
C-CRA             CALL RFFTMLT(U,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
             CALL FFT99M (U,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
C-CRA             CALL RFFTMLT(V,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
             CALL FFT99M (V,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,1)
       RETURN
       END
