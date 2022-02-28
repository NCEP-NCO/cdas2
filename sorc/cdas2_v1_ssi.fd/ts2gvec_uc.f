       SUBROUTINE TS2GVEC(ZS,DS,U,V,JCAP,NLON,NLATH,QLN,RLN,TRIGS,IFAX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    TRANSPOSE OF S2GVEC  
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
C-CRA             DIMENSION WGTS(2*(JCAP+1))
 
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
             DIMENSION WGTS(2*(62+1))
C--------
C-------- INTERNAL SCRATCH DYNAMIC SPACE FOLLOWS:
C--------
C--------
C-CRA                   WGTS=2.*NLON
C          DIMENSION WGTS(2*(JCAP+1))
          DO ITMP=1,2*(JCAP+1)
          WGTS(ITMP)=2.*NLON
          ENDDO
         WGTS(1)=NLON
         WGTS(2)=0.
C-CRA                   DS=0.
C          DIMENSION DS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          DS(ITMP)=0.
          ENDDO
C-CRA                   ZS=0.
C          DIMENSION ZS((JCAP+1)*(JCAP+2))
          DO ITMP=1,(JCAP+1)*(JCAP+2)
          ZS(ITMP)=0.
          ENDDO
C--------
C-------- FIRST DO FOURIER ANALYSIS IN LONGITUDE
C--------
         LOT=NLATH*2
         NLAX=LOT+1
C-CRA             CALL RFFTMLT(U,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
             CALL FFT99M (U,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
C-CRA             CALL RFFTMLT(V,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
             CALL FFT99M (V,WORK,TRIGS,IFAX,NLAX,1,NLON,LOT,-1)
         DO J=1,NLATH
          JR=2*NLATH+1-J
C---------- SEPARATE EVEN AND ODD PARTS
          DO LL=1,2*JCAP+2
           UE(LL)=(U(J,LL)+U(JR,LL))*WGTS(LL)
           UO(LL)=(U(J,LL)-U(JR,LL))*WGTS(LL)
           VE(LL)=(V(J,LL)+V(JR,LL))*WGTS(LL)
           VO(LL)=(V(J,LL)-V(JR,LL))*WGTS(LL)
          END DO
          II0=0
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            ZS(II0+LL)=ZS(II0+LL)+QLN(II0+LL,J)*VE(LL)
            DS(II0+LL)=DS(II0+LL)+QLN(II0+LL,J)*UE(LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             ZS(II0+LL)=ZS(II0+LL)+QLN(II0+LL,J)*VO(LL)
             DS(II0+LL)=DS(II0+LL)+QLN(II0+LL,J)*UO(LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
         END DO
C---------------
C-------------TMULTIPLY VORT, DIV BY I
         DO I=1,(JCAP+1)*(JCAP+2),2
          VORR=ZS(I)
          VORI=ZS(I+1)
          ZS(I)=-VORI
          ZS(I+1)=VORR
          DIVR=DS(I)
          DIVI=DS(I+1)
          DS(I)=-DIVI
          DS(I+1)=DIVR
         END DO
C---------------
         DO J=1,NLATH
          JR=2*NLATH+1-J
C---------- SEPARATE EVEN AND ODD PARTS
          DO LL=1,2*JCAP+2
           UE(LL)=(U(J,LL)+U(JR,LL))*WGTS(LL)
           UO(LL)=(U(J,LL)-U(JR,LL))*WGTS(LL)
           VE(LL)=(V(J,LL)+V(JR,LL))*WGTS(LL)
           VO(LL)=(V(J,LL)-V(JR,LL))*WGTS(LL)
          END DO
          II0=0
          DO M=0,JCAP,2
           DO LL=1,2*(JCAP+1-M)
            ZS(II0+LL)=ZS(II0+LL)+RLN(II0+LL,J)*UO(LL)
            DS(II0+LL)=DS(II0+LL)-RLN(II0+LL,J)*VO(LL)
           END DO
           IF(M.LT.JCAP) THEN
            II0=II0+2*(JCAP+1-M)
            DO LL=1,2*(JCAP-M)
             ZS(II0+LL)=ZS(II0+LL)+RLN(II0+LL,J)*UE(LL)
             DS(II0+LL)=DS(II0+LL)-RLN(II0+LL,J)*VE(LL)
            END DO
            II0=II0+2*(JCAP-M)
           END IF
          END DO
         END DO
       RETURN
       END
