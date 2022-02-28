      FUNCTION ISRCHFGE(N,A,IS,VAL)
      DIMENSION A(N)
      ISRCHFGE=N
      DO NN=IS,N
      IF( A(NN).GE.VAL ) THEN
        ISRCHFGE=NN
        RETURN
      ENDIF
      ENDDO
      RETURN
      END
