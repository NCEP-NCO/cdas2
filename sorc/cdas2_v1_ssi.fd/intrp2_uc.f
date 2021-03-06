      SUBROUTINE INTRP2(F,G,DX,DY,NX,NY,N)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    INTRP2      LINEAR INTERPOLATION IN TWO DIMENSIONS.
C   PRGMMR: PARRISH          ORG: W/NMC22    DATE: 90-10-11
C
C ABSTRACT: LINEAR INTERPOLATE IN 2 DIMS (2ND DIM ALWAYS PERIODIC).
C
C PROGRAM HISTORY LOG:
C   90-10-11  PARRISH
C
C   INPUT ARGUMENT LIST:
C     F        - INPUT INTERPOLATOR
C     DX,DY    - INPUT X,Y -COORDS OF INTERPOLATION POINTS (GRID UNITS)
C     NX,NY    - X,Y-DIMENSIONS OF INTERPOLATOR GRID
C     NY       - DITTO Y
C     N        - NUMBER OF INTERPOLATEES
C
C   OUTPUT ARGUMENT LIST:
C     G        - OUTPUT INTERPOLATEES
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C--------
      DIMENSION F(NX+1,NY+2),G(N),DX(N),DY(N)
C--------
      DO 100 I=1,N
        IX=DX(I)
        IY=DY(I)
        IX=MAX(1,MIN(IX,NX))
        IXP=IX+1
        IXP=MIN(IXP,NX)
        DELX=DX(I)-IX
        DELY=DY(I)-IY
        IF(IY.LT.1) IY=IY+NY
        IF(IY.GT.NY) IY=IY-NY
        IYP=IY+1
        IF(IYP.GT.NY) IYP=IYP-NY
        DELX=MAX(0.,MIN(DELX,1.))
        G(I)=F(IX,IY)*(1.-DELX)*(1.-DELY)
     *      +F(IXP,IY)*DELX*(1.-DELY)
     *      +F(IX,IYP)*(1.-DELX)*DELY
     *      +F(IXP,IYP)*DELX*DELY
100   CONTINUE
      RETURN
      END
