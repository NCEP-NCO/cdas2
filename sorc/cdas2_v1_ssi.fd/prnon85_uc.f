      SUBROUTINE PRNON85(ON85)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    PRNON85   PRINT ON85 INFORMATION
C   PRGMMR: PARRISH        ORG: W/NMC22    DATE: 90-10-10
C
C ABSTRACT: PRINT ON85 INFORMATION.
C
C PROGRAM HISTORY LOG:
C   90-10-10  PARRISH
C
C   INPUT ARGUMENT LIST:
C     ON85     - OFFICE NOTE 85 RECORD                     
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C ATTRIBUTES:
C   LANGUAGE: CFT77
C   MACHINE:  CRAY YMP
C
C$$$
C--------
      INTEGER ON85(4)
      WRITE(6,100)ON85
100   FORMAT(1H ,4Z19)
      RETURN
      END