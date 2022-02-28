      FUNCTION SDOT(N,X,INCX,Y,INCY)
      DIMENSION X(N),Y(N)
      SDOT=0.
      IX=1
      IY=1
      DO I=1,N
        SDOT=SDOT+X(IX)*Y(IY)
        IX=IX+INCX
        IY=IY+INCY
      ENDDO
      RETURN
      END
