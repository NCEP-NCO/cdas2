      SUBROUTINE RDGB(LUGB,LGRIB,LSKIP,KPDS,KGDS,NDATA,LBMS,DATA,LUPTR)
C
C  READ GRIB FILE
C  INPUT
C    LUGB - LOGICAL UNIT TO READ
C    LGRIB - LENGTH OF GRIB RECORD
C    LSKIP - BYTES TO SKIP FOR GRIB RECORD
C  OUTPUT
C    KPDS(25) - UNPACKED PRODUCT DEFINITION SECTION
C    KGDS(22) - UNPACKED GRID DEFINITION SECTION
C    NDATA    - NUMBER OF DATA POINTS
C    LBMS(NDATA) - LOGICAL BIT MAP
C    DATA(NDATA) - DATA UNPACKED
C
      PARAMETER(LLGRIB=720*361)
C     CHARACTER GRIB(LGRIB)*1
      CHARACTER GRIB(LLGRIB)*1
      INTEGER KPDS(25),KGDS(22),KPTR(16)
      LOGICAL LBMS(*)
      REAL DATA(*)
      NDATA=0
        WRITE(*,*) ' rdgb >> baread lugrb=', lugb
      CALL BAREAD(LUGB,LSKIP,LGRIB,LREAD,GRIB)
        WRITE(*,*) ' rdgb <<baread lread,lgrib=',lread, lgrib
      IF(LREAD.LT.LGRIB) THEN
        WRITE(*,*) ' ERROR IN RDGB.  LREAD.LT.LGRIB'
        WRITE(LUPTR,*) ' ERROR IN RDGB.  LREAD.LT.LGRIB'
        CALL ABORT
      ENDIF
      CALL W3FI63(GRIB,KPDS,KGDS,LBMS,DATA,KPTR,IRET)
      IF(IRET.NE.0) THEN
        WRITE(*,*) ' ERROR IN RDGB.  IRET.NE.0 from W3FI63'
        WRITE(LUPTR,*) ' ERROR IN RDGB.  IRET.NE.0 from W3FI63'
        WRITE(LUPTR,*) ' IRET=',IRET
        CALL ABORT
      ENDIF
      NDATA=KPTR(10)
        write(*,*) 'return rdgrb ndata=',ndata
      RETURN
      END
