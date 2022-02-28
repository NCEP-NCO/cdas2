      SUBROUTINE INITW(RU,RV,NLATH,NLON,NSIG,NWRECS,UVFILE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    INITW       SETUP INITIAL RHS FOR WINDS         
C   PRGMMR: DERBER          ORG: W/NMC23    DATE: 91-02-26
C
C ABSTRACT: SETUP INITIAL RHS FOR WINDS
C
C PROGRAM HISTORY LOG:
C   91-02-26  DERBER
C
C   INPUT ARGUMENT LIST:
C     NLATH    - HALF THE NUMBER OF LATITUDES ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NSIG     - NUMBER OF SIGMA LEVELS
C     NWRECS   - NUMBER OF WIND RECORDS
C     NBLK     - BLOCKING FACTOR FOR IUNIT
C     IUNIT    - DATA SCRATCH FILE
C
C   OUTPUT ARGUMENT LIST:
C     RU       - RESULTS FROM OBSERVATION OPERATOR (0 FOR NO DATA)
C     RV       - RESULTS FROM OBSERVATION OPERATOR (0FOR NO DATA)
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION RU(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION RV(2*NLATH+1,NLON+2,NSIG)
C-CRA          DIMENSION UVFILE(*)
 
          DIMENSION RU(2*48+1,192+2,28)
          DIMENSION RV(2*48+1,192+2,28)
          DIMENSION UVFILE(*)
C--------
C-CRA                RU=0.
C       DIMENSION RU(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RU(ITMP,1,1)=0.
          ENDDO
C-CRA                RV=0.
C       DIMENSION RV(2*NLATH+1,NLON+2,NSIG)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)*NSIG
          RV(ITMP,1,1)=0.
          ENDDO
      IF(NWRECS .EQ. 0)RETURN
      NPP=17
C--------
C
C--------
      IS=1
      DO 100 I=1,NWRECS
        NGRP=UVFILE(IS)
        IS=IS+1
      DO 101 K=1,NGRP
        JLAT=UVFILE((K-1)*NPP+IS)
        JLON=UVFILE((K-1)*NPP+IS+1)
        JSIG=UVFILE((K-1)*NPP+IS+2)
        JLATP=UVFILE((K-1)*NPP+IS+3)
        JLONP=UVFILE((K-1)*NPP+IS+4)
        JSIGP=UVFILE((K-1)*NPP+IS+5)
        WGT000=UVFILE((K-1)*NPP+IS+6)
        WGT100=UVFILE((K-1)*NPP+IS+7)
        WGT010=UVFILE((K-1)*NPP+IS+8)
        WGT110=UVFILE((K-1)*NPP+IS+9)
        AERR=UVFILE((K-1)*NPP+IS+16)
        VALU=-UVFILE((K-1)*NPP+IS+14)*AERR
        VALV=-UVFILE((K-1)*NPP+IS+15)*AERR
C       UGES=UVFILE((K-1)*NPP+IS+19)
C       VGES=UVFILE((K-1)*NPP+IS+20)
C       ITYP=UVFILE((K-1)*NPP+IS+21)
C       IF(ITYP .NE. 283)THEN
        WGT001=UVFILE((K-1)*NPP+IS+10)
        WGT101=UVFILE((K-1)*NPP+IS+11)
        WGT011=UVFILE((K-1)*NPP+IS+12)
        WGT111=UVFILE((K-1)*NPP+IS+13)
        RU(JLAT,JLON,JSIG)=RU(JLAT,JLON,JSIG)+WGT000*VALU
        RU(JLATP,JLON,JSIG)=RU(JLATP,JLON,JSIG)+WGT100*VALU
        RU(JLAT,JLONP,JSIG)=RU(JLAT,JLONP,JSIG)+WGT010*VALU
        RU(JLATP,JLONP,JSIG)=RU(JLATP,JLONP,JSIG)+WGT110*VALU
        RU(JLAT,JLON,JSIGP)=RU(JLAT,JLON,JSIGP)+WGT001*VALU
        RU(JLATP,JLON,JSIGP)=RU(JLATP,JLON,JSIGP)+WGT101*VALU
        RU(JLAT,JLONP,JSIGP)=RU(JLAT,JLONP,JSIGP)+WGT011*VALU
        RU(JLATP,JLONP,JSIGP)=RU(JLATP,JLONP,JSIGP)+WGT111*VALU
        RV(JLAT,JLON,JSIG)=RV(JLAT,JLON,JSIG)+WGT000*VALV
        RV(JLATP,JLON,JSIG)=RV(JLATP,JLON,JSIG)+WGT100*VALV
        RV(JLAT,JLONP,JSIG)=RV(JLAT,JLONP,JSIG)+WGT010*VALV
        RV(JLATP,JLONP,JSIG)=RV(JLATP,JLONP,JSIG)+WGT110*VALV
        RV(JLAT,JLON,JSIGP)=RV(JLAT,JLON,JSIGP)+WGT001*VALV
        RV(JLATP,JLON,JSIGP)=RV(JLATP,JLON,JSIGP)+WGT101*VALV
        RV(JLAT,JLONP,JSIGP)=RV(JLAT,JLONP,JSIGP)+WGT011*VALV
        RV(JLATP,JLONP,JSIGP)=RV(JLATP,JLONP,JSIGP)+WGT111*VALV
C       ELSE
C       VALU=WGT000*SU(JLAT,JLON,JSIG)+WGT100*SU(JLATP,JLON,JSIG)
C    *   +WGT010*SU(JLAT,JLONP,JSIG)+WGT110*SU(JLATP,JLONP,JSIG)
C       VALV=WGT000*SV(JLAT,JLON,JSIG)+WGT100*SV(JLATP,JLON,JSIG)
C    *   +WGT010*SV(JLAT,JLONP,JSIG)+WGT110*SV(JLATP,JLONP,JSIG)
C       UANL=UGES*AERR+VALU
C       VANL=VGES*AERR+VALV
C       SPDANL=SQRT(UANL*UANL+VANL*VANL)
C       SPDN=(SPDANL-AERR*UDAT)/SPDANL
C       VALU=AERR*UANL*SPDN
C       VALV=AERR*VANL*SPDN
C       RU(JLAT,JLON,JSIG)=RU(JLAT,JLON,JSIG)+WGT000*VALU
C       RU(JLATP,JLON,JSIG)=RU(JLATP,JLON,JSIG)+WGT100*VALU
C       RU(JLAT,JLONP,JSIG)=RU(JLAT,JLONP,JSIG)+WGT010*VALU
C       RU(JLATP,JLONP,JSIG)=RU(JLATP,JLONP,JSIG)+WGT110*VALU
C       RV(JLAT,JLON,JSIG)=RV(JLAT,JLON,JSIG)+WGT000*VALV
C       RV(JLATP,JLON,JSIG)=RV(JLATP,JLON,JSIG)+WGT100*VALV
C       RV(JLAT,JLONP,JSIG)=RV(JLAT,JLONP,JSIG)+WGT010*VALV
C       RV(JLATP,JLONP,JSIG)=RV(JLATP,JLONP,JSIG)+WGT110*VALV
C       END IF 
 101   CONTINUE
       IS=IS+NGRP*NPP
100   CONTINUE
      RETURN
      END