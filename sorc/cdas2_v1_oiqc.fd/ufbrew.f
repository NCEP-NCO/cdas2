C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE UFBREW(LUNIT,MSG)
 
      CHARACTER*20  MSTR
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      MSG = 0
 
C  MAKE SURE THE FILE IS CLOSED TO COUNT MESSAGES
C  ----------------------------------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.NE.0) GOTO 900
 
C  REWIND AND COUNT THE DATA MESSAGES
C  ----------------------------------
 
      REWIND LUNIT
 
1     READ(LUNIT,END=100,ERR=901) MSTR
      IF(MSTR(1:4).NE.'BUFR'.and.MSTR(1:4).ne.'bufr') GOTO 901
      IF(ICHAR(MSTR(17:17)).NE.11) MSG = MSG+1
      GOTO 1
 
100   CLOSE(LUNIT)
      RETURN
 
C  ERROR EXITS
C  -----------
 
900   WRITE(*,*) 'UFBREW - FILE ALREADY OPEN        '
      STOP 8
901   WRITE(*,*) 'UFBREW - FILE HAS NON-BUFR DATA   '
      STOP 8
      END
