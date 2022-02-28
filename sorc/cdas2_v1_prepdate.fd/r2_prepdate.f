      PROGRAM PREPDATE
 
      CHARACTER*8  SUBSET
 
      DATA LUNIN /11    /
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      iy=0
      im=0
      id=0
      ih=0
      CALL iDATEBF(LUNIN,IY,IM,ID,IH,IDATE)
	  PRINT ('(i4.4,3i2.2)'),IY,IM,ID,IH
	  end

      SUBROUTINE iDATEBF(LUNIN,iy,im,id,ih,IDATE)
	  character*10 CDATE


C  OPEN THE INPUT AND OUTPUT FILES
C  -------------------------------
 
      CALL OPENBF(LUNIN,'IN ',LUNIN)
	  CALL DATELEN(10)
C  -------------------------------------------------------
      if (IREADMG(LUNIN,SUBSET,IDATE).EQ.0) then
			write(cdate,'(I10)') IDATE
			jdate=0
			read (cdate(1:4), '(i4)') iy
			read (cdate(5:6), '(i2)') im
			read (cdate(7:8), '(i2)') id
			read (cdate(9:10),'(i2)') ih
      ENDIF
	  CALL CLOSBF(LUNIN)
 
C  END OF ALL PROCESSING
C  ---------------------
 
      END
