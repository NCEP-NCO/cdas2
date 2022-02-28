      SUBROUTINE SPRUV(WDATA,UGES,VGES,FACT,WTYPE,NWDTA,NWRECS,
     *  NLAT,NLON,NSIG,UVFILE,ERMAX,ERMIN,GROSS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    SPRUV       SAVE WIND INFOR.        
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-13
C
C ABSTRACT: SAVE WIND INFORMATION FOR LATER USE.
C
C PROGRAM HISTORY LOG:
C   90-10-13  PARRISH
C
C   INPUT ARGUMENT LIST:
C     WDATA    - OBS INFO AT OBS LOCATIONS
C     UGES,VGES - GUESS U,V VALUES
C     FACT     - REDUCTION TO 10,20 M FACTOR
C     WTYPE    - WIND TYPE
C     NWDTA    - NUMBER OF OBS
C     NLAT     - NUMBER OF LATITUDES ON GAUSSIAN GRID
C     NLON     - NUMBER OF LONGITUDES ON GAUSSIAN GRID
C     NSIG     - NUMBER OF LAYERS ON GAUSSIAN GRID
C     NBLK     - BLOCKING FACTOR FOR OUTPUT FILE
C     IUNIT    - OUTPUT SCRATCH UNIT
C     ERMAX,ERMIN,GROSS - PARAMETERS FOR GROSS ERROR TEST
C
C   OUTPUT ARGUMENT LIST:
C     NWRECS   - NUMBER OF RECORDS OF WIND DATA
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C
C-CRA          DIMENSION WDATA(NWDTA,7)
C-CRA          DIMENSION UGES(NWDTA),VGES(NWDTA),WTYPE(NWDTA)
C-CRA          DIMENSION NUMW(NSIG)
C-CRA          DIMENSION FACT(NWDTA)
C-CRA          DIMENSION UPLTY(NSIG),VPLTY(NSIG)
C-CRA          DIMENSION ICOUNT(NLAT,NLON,NSIG),NCOUNT(NWDTA)
 
C          DIMENSION WDATA(_NWDTA_,7)
          DIMENSION WDATA(85000,7)
C          DIMENSION UGES(_NWDTA_),VGES(_NWDTA_),WTYPE(_NWDTA_)
          DIMENSION UGES(85000),VGES(85000),WTYPE(85000)
          DIMENSION NUMW(28)
C          DIMENSION FACT(_NWDTA_)
          DIMENSION FACT(85000)
          DIMENSION UPLTY(28),VPLTY(28)
C          DIMENSION ICOUNT(96,192,28),NCOUNT(_NWDTA_)
          DIMENSION ICOUNT(96,192,28),NCOUNT(85000)

           DIMENSION JL(128,6)
           DIMENSION UVFILE(*)
C--------
      OBERMAX=-1.E50
      OBERMIN=1.E50
      RESMAX=-1.E50
      RATMAX=-1.E50
      NUMGROSS=0
C-CRA                NUMW=0
C       DIMENSION NUMW(NSIG)
          DO ITMP=1,NSIG
          NUMW(ITMP)=0
          ENDDO
      NTOT=0
      IER=1
      ILON=2
      ILAT=3
      ISIG=4
      IURES=5
      IVRES=6
      NLTH=NLAT/2
      FACTOR=2.0
      SSMPEN=0.
      NUMSSM=0
C-CRA                UPLTY=0.
C       DIMENSION UPLTY(NSIG),VPLTY(NSIG)
          DO ITMP=1,NSIG
          UPLTY(ITMP)=0.
          ENDDO
      UMPLTY=0.
C-CRA                VPLTY=0.
C       DIMENSION UPLTY(NSIG),VPLTY(NSIG)
          DO ITMP=1,NSIG
          VPLTY(ITMP)=0.
          ENDDO
      VMPLTY=0.
      NPP=17
      ANLON=FLOAT(NLON)
      IREC=0
      INC=NWDTA/128
      INC=MAX(INC,1)
      NWTTOT=0
C-CRA                NCOUNT=0
C       DIMENSION ICOUNT(NLAT,NLON,NSIG),NCOUNT(NWDTA)
          DO ITMP=1,NWDTA
          NCOUNT(ITMP)=0
          ENDDO
      IS=1
      DO 200 KK=1,NWDTA
        NGRP=0
C-CRA                  ICOUNT=0
C       DIMENSION ICOUNT(NLAT,NLON,NSIG),NCOUNT(NWDTA)
          DO ITMP=1,NLAT*NLON*NSIG
          ICOUNT(ITMP,1,1)=0
          ENDDO
        ISSAVE=IS
        IS=IS+1
        I128=1
        NUMDAT=0
        DO 100 III=1,5*INC
        IBEG=MOD(III-1,INC)+1
        DO 100 I=IBEG,NWDTA,INC
        IF(NCOUNT(I) .GT. 0)GO TO 100
        JLAT=WDATA(I,ILAT)
        IF(WDATA(I,ILON).GE. ANLON+1.)WDATA(I,ILON)=
     *     WDATA(I,ILON)-ANLON
        IF(WDATA(I,ILON).LT. 1.)WDATA(I,ILON)=WDATA(I,ILON)+ANLON
        JLON=WDATA(I,ILON)
        JSIG=WDATA(I,ISIG)
        DX=WDATA(I,ILON)-JLON
        DY=WDATA(I,ILAT)-JLAT
        DS=WDATA(I,ISIG)-JSIG
        JLAT=MAX(1,MIN(JLAT,NLAT))
        JSIG=MAX(1,MIN(JSIG,NSIG))
        IF(ICOUNT(JLAT,JLON,JSIG) .EQ. 1)GO TO 100
        JLATP=JLAT+1
        JLATP=MIN(JLATP,NLAT)
        IF(ICOUNT(JLATP,JLON,JSIG) .EQ. 1)GO TO 100
        JSIGP=JSIG+1
        JSIGP=MIN(JSIGP,NSIG)
        IF(ICOUNT(JLAT,JLON,JSIGP) .EQ. 1)GO TO 100
        IF(ICOUNT(JLATP,JLON,JSIG) .EQ. 1)GO TO 100
        JLONP=JLON+1
        IF(JLONP .GT. NLON)JLONP=JLONP-NLON
        IF(ICOUNT(JLAT,JLONP,JSIGP) .EQ. 1)GO TO 100
        IF(ICOUNT(JLAT,JLONP,JSIG) .EQ. 1)GO TO 100
        IF(ICOUNT(JLATP,JLONP,JSIGP) .EQ. 1)GO TO 100
        IF(ICOUNT(JLATP,JLONP,JSIG) .EQ. 1)GO TO 100
        IF(JLAT .LE. NLTH)WDATA(I,IER)=WDATA(I,IER)*FACTOR
C-----------------------------------GROSS ERROR TEST ADDED HERE
        OBSERROR=1./MAX(SQRT(WDATA(I,IER)),1.E-10)
        OBSERRLM=MAX(ERMIN,MIN(ERMAX,OBSERROR))
        RESIDUAL=SQRT(WDATA(I,IURES)**2+WDATA(I,IVRES)**2)
        RATIO=RESIDUAL/OBSERRLM
        IF(OBSERROR.LT.1.E5) OBERMAX=MAX(OBERMAX,OBSERROR)
        OBERMIN=MIN(OBERMIN,OBSERROR)
        RESMAX=MAX(RESMAX,RESIDUAL)
        RATMAX=MAX(RATMAX,RATIO)
        IF(RATIO.GT.GROSS) THEN
         NUMGROSS=NUMGROSS+1
         WDATA(I,IER)=0.
        END IF
        IF(NINT(WTYPE(I)) .NE. 283)THEN
        UPLTY(JSIG)=UPLTY(JSIG)+WDATA(I,IER)*WDATA(I,IURES)**2
        VPLTY(JSIG)=VPLTY(JSIG)+WDATA(I,IER)*WDATA(I,IVRES)**2
        NUMW(JSIG)=NUMW(JSIG)+1
        ELSE
        SSMPEN=SSMPEN+WDATA(I,IER)*(WDATA(I,IURES)-WDATA(I,IVRES))
     *     **2
        NUMSSM=NUMSSM+1
        END IF
        ICOUNT(JLAT,JLON,JSIG)=1
        ICOUNT(JLATP,JLON,JSIG)=1
        ICOUNT(JLAT,JLONP,JSIG)=1
        ICOUNT(JLATP,JLONP,JSIG)=1
        ICOUNT(JLAT,JLON,JSIGP)=1
        ICOUNT(JLATP,JLON,JSIGP)=1
        ICOUNT(JLAT,JLONP,JSIGP)=1
        ICOUNT(JLATP,JLONP,JSIGP)=1
        JL(I128,1)=JLAT
        JL(I128,2)=JLON
        JL(I128,3)=JSIG
        JL(I128,4)=JLATP
        JL(I128,5)=JLONP
        JL(I128,6)=JSIGP
        WDATA(I,IER)=SQRT(WDATA(I,IER))
        WGT000=FACT(I)*WDATA(I,IER)*(1.0-DX)*(1.0-DY)*(1.0-DS)
        WGT010=FACT(I)*WDATA(I,IER)*DX*(1.0-DY)*(1.0-DS)
        WGT100=FACT(I)*WDATA(I,IER)*(1.-DX)*DY*(1.0-DS)
        WGT110=FACT(I)*WDATA(I,IER)*DX*DY*(1.0-DS)
        WGT001=FACT(I)*WDATA(I,IER)*(1.-DX)*(1.-DY)*DS
        WGT011=FACT(I)*WDATA(I,IER)*DX*(1.-DY)*DS
        WGT101=FACT(I)*WDATA(I,IER)*(1.-DX)*DY*DS
        WGT111=FACT(I)*WDATA(I,IER)*DX*DY*DS
        UVFILE(IS)=JLAT+.001
        UVFILE(IS+1)=JLON+.001
        UVFILE(IS+2)=JSIG+.001
        UVFILE(IS+3)=JLATP+.001
        UVFILE(IS+4)=JLONP+.001
        UVFILE(IS+5)=JSIGP+.001
        UVFILE(IS+6)=WGT000
        UVFILE(IS+7)=WGT100
        UVFILE(IS+8)=WGT010
        UVFILE(IS+9)=WGT110
        UVFILE(IS+10)=WGT001
        UVFILE(IS+11)=WGT101
        UVFILE(IS+12)=WGT011
        UVFILE(IS+13)=WGT111
        UVFILE(IS+14)=WDATA(I,IURES)
        UVFILE(IS+15)=WDATA(I,IVRES)
        UVFILE(IS+16)=WDATA(I,IER)
C       UVFILE(IS+17)=UGES(I)
C       UVFILE(IS+18)=VGES(I)
C       UVFILE(IS+19)=WTYPE(I)+.001
        NWTTOT=NWTTOT+1
        IS=IS+NPP
        NGRP=NGRP+1
        I128=I128+1
        NCOUNT(I)=1
        IF(I128 .EQ. 129)I128=1
        IF(NGRP .GT. 128)THEN
        ICOUNT(JL(I128,1),JL(I128,2),JL(I128,3))=0
        ICOUNT(JL(I128,4),JL(I128,2),JL(I128,3))=0
        ICOUNT(JL(I128,1),JL(I128,5),JL(I128,3))=0
        ICOUNT(JL(I128,4),JL(I128,5),JL(I128,3))=0
        ICOUNT(JL(I128,1),JL(I128,2),JL(I128,6))=0
        ICOUNT(JL(I128,4),JL(I128,2),JL(I128,6))=0
        ICOUNT(JL(I128,1),JL(I128,5),JL(I128,6))=0
        ICOUNT(JL(I128,4),JL(I128,5),JL(I128,6))=0
        END IF
100   CONTINUE
121   CONTINUE
      UVFILE(ISSAVE)=NGRP+.001
      IREC=IREC+1
      IF(NWTTOT .EQ. NWDTA)GO TO 201
 200  CONTINUE
 201  NWRECS=IREC
      DO 251 K=1,NSIG
        UMPLTY=UMPLTY+UPLTY(K)
        VMPLTY=VMPLTY+VPLTY(K)
        NTOT=NTOT+NUMW(K)
        WRITE(6,240)NUMW(K),K,UPLTY(K),VPLTY(K)
240     FORMAT(' NUMBER OF WIND OBS=',I9,' AT LEVEL ',I4,' PEN = ',
     *     2E12.4)
251   CONTINUE
      PRINT *,' TOTAL NUMBER OF SSM/I OBS=',NUMSSM
      PRINT *,' TOTAL NUMBER OF WIND OBS=',NTOT
      PRINT *,' SSM/I PENALTY = ', SSMPEN
      TU=UMPLTY/FLOAT(NTOT)
      WRITE(6,930)UMPLTY,TU
930   FORMAT(' TOTAL U OBS PENALTY=',E12.4,E12.4)
      TV=VMPLTY/FLOAT(NTOT)
      WRITE(6,940)VMPLTY,TV
940   FORMAT(' TOTAL V OBS PENALTY=',E12.4,E12.4)
       WRITE(6,*)' GROSS ERROR CHECK FOR WINDS:'
       WRITE(6,*)'   OBS ERROR MAX,MIN=',OBERMAX,OBERMIN
       WRITE(6,*)'   FOR CHECK, OBS ERROR BOUNDED BY ',ERMIN,ERMAX
       WRITE(6,*)'   FOR CHECK, MAX RATIO RESIDUAL/OB ERROR =',GROSS
       WRITE(6,*)'   MAX RESIDUAL=',RESMAX
       WRITE(6,*)'   MAX RATIO=',RATMAX
       WRITE(6,*)'   NUMBER OBS THAT FAILED GROSS TEST = ',NUMGROSS
      RETURN
      END