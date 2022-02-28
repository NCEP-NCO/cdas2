      SUBROUTINE INTQPW(RT,NLATH,NLON,NSIG,NQRECS,NPWRECS,
     *     PWCON,QFILE,PWFILE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    INTQPW       APPLY OBSERVATION OPERATOR FOR Q AND P.W.
C   PRGMMR: DERBER          ORG: W/NMC23    DATE: 91-02-26
C
C ABSTRACT: APPLY OBSERVATION OPERATOR FOR Q AND PRECIP. WATER
C
C PROGRAM HISTORY LOG:
C   91-02-26  DERBER
C
C   INPUT ARGUMENT LIST:
C     RT       - SEARCH DIRECTION FOR Q
C     NLATH    - HALF THE NUMBER OF LATITUDES ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NSIG     - NUMBER OF SIGMA LEVELS
C     NQRECS   - NUMBER OF Q RECORDS
C     NPWRECS  - NUMBER OF PRECIP. WATER RECORDS
C     PWCON    - VERTICAL INTEGRATION PRECIP. WATER CONSTANTS
C     NBLK     - BLOCKING FACTOR FOR IUNIT
C     IUNIT    - DATA SCRATCH FILE
C
C   OUTPUT ARGUMENT LIST:
C     RT       - RESULTS FROM OBSERVATION OPERATOR (0 FOR NO DATA)
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION PWCON(NSIG)
C-CRA          DIMENSION RT(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION QFILE(*),PWFILE(*)
C-CRA          DIMENSION ST(2*NLATH+1,NLON+2,NSIG)
 
          DIMENSION PWCON(28)
          DIMENSION RT(2*48+1,192+2,28)
          DIMENSION QFILE(*),PWFILE(*)
          DIMENSION ST(2*48+1,192+2,28)
C--------
C--------
C-------- INITIALIZE GRIDS
C--------
C-CRA                ST=RT
C       DIMENSION ST(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          ST(ITMP,1,1)=RT(ITMP,1,1)
          ENDDO
C-CRA                RT=0.
C       DIMENSION RT(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RT(ITMP,1,1)=0.
          ENDDO
      IF(NPWRECS .EQ. 0)GO TO 1000
      NPP=11
      IS=1
      DO 700 I=1,NPWRECS
        NGRP=PWFILE(IS)
        IS=IS+1
        DO 701 IRPT=1,NGRP
        JLAT=PWFILE((IRPT-1)*NPP+IS)
        JLON=PWFILE((IRPT-1)*NPP+IS+1)
        JLATP=PWFILE((IRPT-1)*NPP+IS+2)
        JLONP=PWFILE((IRPT-1)*NPP+IS+3)
        WGT00=PWFILE((IRPT-1)*NPP+IS+4)
        WGT10=PWFILE((IRPT-1)*NPP+IS+5)
        WGT01=PWFILE((IRPT-1)*NPP+IS+6)
        WGT11=PWFILE((IRPT-1)*NPP+IS+7)
        PSFC=PWFILE((IRPT-1)*NPP+IS+8)
C       PDAT=PWFILE((IRPT-1)*NPP+IS+9)
C       AERR=PWFILE((IRPT-1)*NPP+IS+10)
C       PWGE=PWFILE((IRPT-1)*NPP+IS+11)
C       PWTY=PWFILE((IRPT-1)*NPP+IS+12)
        VAL=0.
        DO 400 K=1,NSIG
        VAL=(WGT00*ST(JLAT,JLON,K)+WGT10*ST(JLATP,JLON,K)
     *     +WGT01*ST(JLAT,JLONP,K)+WGT11*ST(JLATP,JLONP,K))*
     *     PWCON(K)+VAL
 400    CONTINUE
        VAL=VAL*PSFC*PSFC
C       VAL=(VAL-PDAT*AERR)*PSFC*PSFC
        DO 401 K=1,NSIG
        RT(JLAT,JLON,K)=RT(JLAT,JLON,K)+WGT00*VAL*PWCON(K)
        RT(JLATP,JLON,K)=RT(JLATP,JLON,K)+WGT10*VAL*PWCON(K)
        RT(JLAT,JLONP,K)=RT(JLAT,JLONP,K)+WGT01*VAL*PWCON(K)
        RT(JLATP,JLONP,K)=RT(JLATP,JLONP,K)+WGT11*VAL*PWCON(K)
 401    CONTINUE
 701  CONTINUE
      IS=IS+NPP*NGRP
700   CONTINUE
1000  IF(NQRECS .EQ. 0)RETURN
      IS=1
      NPP=16
C--------
C
C--------
      DO 100 I=1,NQRECS
        NGRP=QFILE(IS)
        IS=IS+1
      DO 101 K=1,NGRP
        JLAT=QFILE((K-1)*NPP+IS)
        JLON=QFILE((K-1)*NPP+IS+1)
        JSIG=QFILE((K-1)*NPP+IS+2)
        JLATP=QFILE((K-1)*NPP+IS+3)
        JLONP=QFILE((K-1)*NPP+IS+4)
        JSIGP=QFILE((K-1)*NPP+IS+5)
        WGT000=QFILE((K-1)*NPP+IS+6)
        WGT100=QFILE((K-1)*NPP+IS+7)
        WGT010=QFILE((K-1)*NPP+IS+8)
        WGT110=QFILE((K-1)*NPP+IS+9)
        WGT001=QFILE((K-1)*NPP+IS+10)
        WGT101=QFILE((K-1)*NPP+IS+11)
        WGT011=QFILE((K-1)*NPP+IS+12)
        WGT111=QFILE((K-1)*NPP+IS+13)
C       TDAT=QFILE((K-1)*NPP+IS+14)
C       AERR=QFILE((K-1)*NPP+IS+15)
C       QGES=QFILE((K-1)*NPP+IS+16)
C       QTYP=QFILE((K-1)*NPP+IS+17)
        VAL=WGT000*ST(JLAT,JLON,JSIG)+WGT100*ST(JLATP,JLON,JSIG)
     *   +WGT010*ST(JLAT,JLONP,JSIG)+WGT110*ST(JLATP,JLONP,JSIG)
     *   +WGT001*ST(JLAT,JLON,JSIGP)+WGT101*ST(JLATP,JLON,JSIGP)
     *   +WGT011*ST(JLAT,JLONP,JSIGP)+WGT111*ST(JLATP,JLONP,JSIGP)
        RT(JLAT,JLON,JSIG)=RT(JLAT,JLON,JSIG)+WGT000*VAL
        RT(JLATP,JLON,JSIG)=RT(JLATP,JLON,JSIG)+WGT100*VAL
        RT(JLAT,JLONP,JSIG)=RT(JLAT,JLONP,JSIG)+WGT010*VAL
        RT(JLATP,JLONP,JSIG)=RT(JLATP,JLONP,JSIG)+WGT110*VAL
        RT(JLAT,JLON,JSIGP)=RT(JLAT,JLON,JSIGP)+WGT001*VAL
        RT(JLATP,JLON,JSIGP)=RT(JLATP,JLON,JSIGP)+WGT101*VAL
        RT(JLAT,JLONP,JSIGP)=RT(JLAT,JLONP,JSIGP)+WGT011*VAL
        RT(JLATP,JLONP,JSIGP)=RT(JLATP,JLONP,JSIGP)+WGT111*VAL
 101   CONTINUE
      IS=IS+NPP*NGRP
100   CONTINUE
      RETURN
      END
