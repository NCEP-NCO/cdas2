      SUBROUTINE PREPS(DLONS,DLATS,DPRES,RSPRES,
     *  MSDAT,PSG,NLAT,NLON,NSIG,GLATS,GLONS,SIGL)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    PREPS       PRELIMINARY STUFF BEFORE RES. CALC. T
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-11
C
C ABSTRACT: PRELIMINARY STUFF BEFORE RESIDUAL CALCULATION FOR SAT TEMPS.
C
C PROGRAM HISTORY LOG:
C   90-10-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     DLONS,DLATS - OBS LONGITUDES AND LATITUDES (RADIANS IN AND OUT)
C     DPRES    - PRES (MB*10+QM IN, GRID COORDS IN SIGMA OUT)
C     MTDAT    - NUMBER OF OBSERVATIONS
C     PSG      - MODEL GUESS LOG(PSFC), P IN CB
C     NLAT     - NUMBER OF GAUSSIAN LATS POLE TO POLE
C     NLON     - NUMBER OF LONGITUDES
C     NSIG     - NUMBER OF SIGMA LEVELS
C     GLATS,GLONS - GRID LATITUDES AND LONGITUDES
C     SIGL     - SIGMA LAYER MIDPOINT VALUES
C
C   OUTPUT ARGUMENT LIST:
C     DLONS,DLATS - OBS LONGITUDINAL AND LATITUDINAL GRID LOCATION 
C     RSPRES   - OBSERVATION PRESSURES
C     AND AS INDICATED ABOVE
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C--------
C
C-CRA          DIMENSION DLONS(MSDAT)
C-CRA          DIMENSION DLATS(MSDAT),DPRES(MSDAT)
C-CRA          DIMENSION PSG(NLAT+1,NLON+2)
C-CRA          DIMENSION GLATS(NLAT),GLONS(NLON),SIGL(NSIG)
C-CRA          DIMENSION RSPRES(MSDAT)
C-CRA          DIMENSION RBPRES(MSDAT)
C-CRA          DIMENSION SIGLL(NSIG)
 
C          DIMENSION DLONS(_MSDAT_)
          DIMENSION DLONS(120000)
C          DIMENSION DLATS(_MSDAT_),DPRES(_MSDAT_)
          DIMENSION DLATS(120000),DPRES(120000)
          DIMENSION PSG(96+1,192+2)
          DIMENSION GLATS(96),GLONS(192),SIGL(28)
C          DIMENSION RSPRES(_MSDAT_)
          DIMENSION RSPRES(120000)
C          DIMENSION RBPRES(_MSDAT_)
          DIMENSION RBPRES(120000)
          DIMENSION SIGLL(28)
C--------
C-------- LOCAL SPACE
C--------
C--------
C-------- GET LOG(SIG)
C--------
C-CRA                SIGLL=LOG(SIGL)
C       DIMENSION SIGLL(NSIG)
          DO ITMP=1,NSIG
          SIGLL(ITMP)=LOG(SIGL(ITMP))
          ENDDO
C--------
C-------- CONVERT OBS LATS AND LONS TO GRID COORDINATES
C--------
      CALL GDCRDP(DLATS,MSDAT,GLATS,NLAT)
      CALL GDCRDP(DLONS,MSDAT,GLONS,NLON)
C--------
C-------- 3.  INTERPOLATE SURFACE PRESSURE
C--------
C-------- OBTAIN GUESS SURFACE PRESSURE AT OBS LOCATIONS
C--------
      CALL INTRP2(PSG,RBPRES,DLATS,DLONS,
     *  NLAT,NLON,MSDAT)
C--------
C-------- CONVERT OBS PRESSURE TO SIGMA, THEN GET GRID COORDINATES
C--------
C-CRA                RSPRES=10.*EXP(DPRES)
C       DIMENSION RSPRES(MSDAT)
          DO ITMP=1,MSDAT
          RSPRES(ITMP)=10.*EXP(DPRES(ITMP))
          ENDDO
C-CRA                DPRES=DPRES-RBPRES
C       DIMENSION DLATS(MSDAT),DPRES(MSDAT)
          DO ITMP=1,MSDAT
          DPRES(ITMP)=DPRES(ITMP)-RBPRES(ITMP)
          ENDDO
      CALL GDCRDN(DPRES,MSDAT,SIGLL,NSIG)
      RETURN
      END
