      SUBROUTINE INTPS(PS,NLATH,NLON,NPRECS,PSFILE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    INTPS       APPLY OBSERVATION OPERATOR FOR PS 
C   PRGMMR: DERBER          ORG: W/NMC23    DATE: 91-02-26
C
C ABSTRACT: APPLY OBSERVATION OPERATOR FOR PS OBSERVATIONS
C
C PROGRAM HISTORY LOG:
C   91-02-26  DERBER
C
C   INPUT ARGUMENT LIST:
C     PS       - SEARCH DIRECTION FOR PS
C     NLATH    - HALF THE NUMBER OF LATITUDES ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NPRECS   - NUMBER OF PS RECORDS
C     NBLK     - BLOCKING FACTOR FOR IUNIT
C     IUNIT    - DATA SCRATCH FILE
C
C   OUTPUT ARGUMENT LIST:
C     PS       - RESULTS FROM OBSERVATION OPERATOR (0 FOR NO DATA)
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION PS(2*NLATH+1,NLON+2)
C-CRA          DIMENSION PSFILE(*)
C-CRA          DIMENSION SPS(2*NLATH+1,NLON+2)
 
          DIMENSION PS(2*48+1,192+2)
          DIMENSION PSFILE(*)
          DIMENSION SPS(2*48+1,192+2)
C--------
C--------
C-------- INITIALIZE GRIDS
C--------
C-CRA                SPS=PS
C       DIMENSION SPS(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          SPS(ITMP,1)=PS(ITMP,1)
          ENDDO
C-CRA                PS=0.
C       DIMENSION PS(2*NLATH+1,NLON+2)
          DO ITMP=1,(2*NLATH+1)*(NLON+2)
          PS(ITMP,1)=0.
          ENDDO
      IF(NPRECS .EQ. 0)RETURN
      NPP=10
      IS=1
      DO 100 I=1,NPRECS
        NGRP=PSFILE(IS)+.001
        IS=IS+1
        DO 101 IRPT=1,NGRP
        JLAT=PSFILE((IRPT-1)*NPP+IS)
        JLON=PSFILE((IRPT-1)*NPP+IS+1)
        JLATP=PSFILE((IRPT-1)*NPP+IS+2)
        JLONP=PSFILE((IRPT-1)*NPP+IS+3)
        WGT00=PSFILE((IRPT-1)*NPP+IS+4)
        WGT10=PSFILE((IRPT-1)*NPP+IS+5)
        WGT01=PSFILE((IRPT-1)*NPP+IS+6)
        WGT11=PSFILE((IRPT-1)*NPP+IS+7)
C       PDAT=PSFILE((IRPT-1)*NPP+IS+8)
C       AERR=PSFILE((IRPT-1)*NPP+IS+9)
C       PGES=PSFILE((IRPT-1)*NPP+IS+10)
C       PTYP=PSFILE((IRPT-1)*NPP+IS+11)
        VAL=WGT00*SPS(JLAT,JLON)+WGT10*SPS(JLATP,JLON)
     *     +WGT01*SPS(JLAT,JLONP)+WGT11*SPS(JLATP,JLONP)
        PS(JLAT,JLON)=PS(JLAT,JLON)+WGT00*VAL
        PS(JLATP,JLON)=PS(JLATP,JLON)+WGT10*VAL
        PS(JLAT,JLONP)=PS(JLAT,JLONP)+WGT01*VAL
        PS(JLATP,JLONP)=PS(JLATP,JLONP)+WGT11*VAL
 101  CONTINUE
      IS=IS+NPP*NGRP
100   CONTINUE
      RETURN
      END
