      SUBROUTINE SPRP(PDATA,PGES,PTYPE,NPSDTA,NPSRECS,NLAT,NLON,
     *   PSFILE,ERMAX,ERMIN,GROSS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    SPRP        STORE PSFC INFO.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-12
C
C ABSTRACT: STORE INFORMATION FOR PSFC DATA.
C
C PROGRAM HISTORY LOG:
C   90-10-12  PARRISH
C
C   INPUT ARGUMENT LIST:
C     PDATA    - OBS INFO AT OBS LOCATIONS
C     PGES     - GUESS VALUES FOR PSFC.
C     PTYPE    - OBSERVATION TYPES
C     NPSDTA   - NUMBER OF OBS
C     NLAT     - NUMBER OF LATITUDES ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NBLK     - BLOCKING DATA FOR FILE IUNIT
C     IUNIT    - OUTPUT FILE FOR OBS. INFORMATION
C     ERMAX,ERMIN,GROSS - PARAMETERS FOR GROSS ERROR CHECK
C
C   OUTPUT ARGUMENT LIST:
C     NPSDTA   - NUMBER OF RECORDS OF PSFC. OBS
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION PDATA(NPSDTA,7),PGES(NPSDTA),PTYPE(NPSDTA)
C-CRA          DIMENSION ICOUNT(NLAT,NLON)
C-CRA          DIMENSION NCOUNT(NPSDTA),JL(128,4)
 
C          DIMENSION PDATA(_NPSDTA_,7),PGES(_NPSDTA_)
          DIMENSION PDATA(18000,7),PGES(18000)
C          DIMENSION PTYPE(_NPSDTA_)
          DIMENSION PTYPE(18000)
          DIMENSION ICOUNT(96,192)
C          DIMENSION NCOUNT(_NPSDTA_),JL(128,4)
          DIMENSION NCOUNT(18000),JL(128,4)

           DIMENSION PSFILE(*)
C--------
      OBERMAX=-1.E50
      OBERMIN=1.E50
      RESMAX=-1.E50
      RATMAX=-1.E50
      NUMGROSS=0
C-CRA                NCOUNT=0
C       DIMENSION NCOUNT(NPSDTA),JL(128,4)
          DO ITMP=1,NPSDTA
          NCOUNT(ITMP)=0
          ENDDO
      IER=1
      ILON=2
      ILAT=3
      IPRES=4
      NSUPERP=0
      PSPLTY=0.
      NUMW=0
      NPP=10
      NDCNT=0
C
C   INCREASE OBS ERRORS IN THE N.H. BY A FACTOR OF 2 TO ACCOUNT
C--------
C-------- INITIALIZE GRIDS
C--------
      ANLON=FLOAT(NLON)
      NLTH=NLAT/2
      FACTOR=2.0
      NRECS=0
      INC=NPSDTA/128
      INC=MAX(INC,1)
      IS=1
      DO 200 K=1,NPSDTA
      ISSAVE=IS
      IS=IS+1
C-CRA                ICOUNT=0
C       DIMENSION ICOUNT(NLAT,NLON)
          DO ITMP=1,NLAT*NLON
          ICOUNT(ITMP,1)=0
          ENDDO
      NUMDAT=0
      I128=1
      DO 100 KK=1,INC*5
      IBEG=MOD(KK-1,INC)+1
      IEND=NPSDTA
      DO 100 I=IBEG,IEND,INC
        IF(NCOUNT(I) .GT. 0)GO TO 100
        JLAT=PDATA(I,ILAT)
        IF(PDATA(I,ILON) .GE. ANLON+1.)PDATA(I,ILON)=
     *        PDATA(I,ILON)-ANLON
        IF(PDATA(I,ILON) .LT. 1.)PDATA(I,ILON)=PDATA(I,ILON)+ANLON
        JLON=PDATA(I,ILON)
        JLAT=MAX(1,MIN(JLAT,NLAT))
        IF(ICOUNT(JLAT,JLON) .GT. 0)GO TO 100
        DY=PDATA(I,ILAT)-JLAT
        DX=PDATA(I,ILON)-JLON
        JLONP=JLON+1
        IF(JLONP .GT. NLON)JLONP=JLONP-NLON
        IF(ICOUNT(JLAT,JLONP) .GT. 0)GO TO 100
        JLATP=JLAT+1
        JLATP=MIN(JLATP,NLAT)
        IF(ICOUNT(JLATP,JLON) .GT. 0)GO TO 100
        IF(ICOUNT(JLATP,JLONP) .GT. 0)GO TO 100
        ICOUNT(JLAT,JLON)=1
        ICOUNT(JLATP,JLON)=1
        ICOUNT(JLAT,JLONP)=1
        ICOUNT(JLATP,JLONP)=1
        JL(I128,1)=JLAT
        JL(I128,2)=JLATP
        JL(I128,3)=JLON
        JL(I128,4)=JLONP
        I128=I128+1
        NCOUNT(I)=1
        NDCNT=NDCNT+1
        IF(JLAT .LE. NLTH)PDATA(I,IER)=PDATA(I,IER)*FACTOR
        PDATA(I,IER)=SQRT(PDATA(I,IER))
C-----------------------------------GROSS ERROR TEST ADDED HERE
        OBSERROR=1000./MAX(PDATA(I,IER),1.E-10)
        OBSERRLM=MAX(ERMIN,MIN(ERMAX,OBSERROR))
        RESIDUAL=ABS(1000.*PDATA(I,IPRES))
        RATIO=RESIDUAL/OBSERRLM
        IF(OBSERROR.LT.1.E5) OBERMAX=MAX(OBERMAX,OBSERROR)
        OBERMIN=MIN(OBERMIN,OBSERROR)
        RESMAX=MAX(RESMAX,RESIDUAL)
        RATMAX=MAX(RATMAX,RATIO)
        IF(RATIO.GT.GROSS) THEN
         NUMGROSS=NUMGROSS+1
         PDATA(I,IER)=0.
        END IF
        VAL=PDATA(I,IER)*PDATA(I,IPRES)
        WGT00=PDATA(I,IER)*(1.0-DX)*(1.0-DY)
        WGT10=PDATA(I,IER)*(1.0-DX)*DY
        WGT01=PDATA(I,IER)*DX*(1.0-DY)
        WGT11=PDATA(I,IER)*DX*DY
        NUMDAT=NUMDAT+1
        PSFILE(IS)=JLAT
        PSFILE(IS+1)=JLON
        PSFILE(IS+2)=JLATP
        PSFILE(IS+3)=JLONP
        PSFILE(IS+4)=WGT00
        PSFILE(IS+5)=WGT10
        PSFILE(IS+6)=WGT01
        PSFILE(IS+7)=WGT11
        PSFILE(IS+8)=PDATA(I,IPRES)
        PSFILE(IS+9)=PDATA(I,IER)
C       PSFILE(IS+10)=PGES(I)
C       PSFILE(IS+11)=PTYPE(I)
        IS=IS+NPP
        NSUPERP=NSUPERP+1
        PSPLTY=PSPLTY+VAL*VAL
        NUMW=NUMW+1
        IF(I128 .EQ. 129)I128=1
        IF(NUMDAT .GT. 128)THEN
          ICOUNT(JL(I128,1),JL(I128,3))=0
          ICOUNT(JL(I128,1),JL(I128,4))=0
          ICOUNT(JL(I128,2),JL(I128,3))=0
          ICOUNT(JL(I128,2),JL(I128,4))=0
        END IF
100   CONTINUE
101   CONTINUE
      PSFILE(ISSAVE)=NUMDAT+.001
      NRECS=NRECS+1
      IF(NDCNT .EQ. NPSDTA)GO TO 201
200   CONTINUE
      PW=PSPLTY/FLOAT(NUMW)
201   WRITE(6,955)PSPLTY,NUMW,PW
955   FORMAT(' TOTAL PSFC OBS PENALTY=',E12.4,I8,E12.4)
      NPSRECS=NRECS
       WRITE(6,*)' GROSS ERROR CHECK FOR PSFCS:'
       WRITE(6,*)'   OBS ERROR MAX,MIN=',OBERMAX,OBERMIN
       WRITE(6,*)'   FOR CHECK, OBS ERROR BOUNDED BY ',ERMIN,ERMAX
       WRITE(6,*)'   FOR CHECK, MAX RATIO RESIDUAL/OB ERROR =',GROSS
       WRITE(6,*)'   MAX RESIDUAL=',RESMAX
       WRITE(6,*)'   MAX RATIO=',RATMAX
       WRITE(6,*)'   NUMBER OBS THAT FAILED GROSS TEST = ',NUMGROSS
      RETURN
      END
