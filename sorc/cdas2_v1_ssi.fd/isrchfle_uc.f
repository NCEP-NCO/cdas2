      FUNCTION ISRCHFLE(N,A,IS,VAL)
      DIMENSION A(N)
      ISRCHFLE=N
      DO NN=IS,N
      IF( A(NN).LE.VAL ) THEN
        ISRCHFLE=NN
        RETURN
      ENDIF
      ENDDO
      RETURN
      END
