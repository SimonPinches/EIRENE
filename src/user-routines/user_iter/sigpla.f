C 
C
      SUBROUTINE EIRENE_SIGPLA(INIT,JJJ,ZDS,DUM1,PSIG,DUM2,ARGST)
C  MG 27/8/21 (MG) User routine to integrate plasma parameters,
C  e.g. electron density    
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C               
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO
C          THE PARAMETER PSIG (TOTAL) and ARGST (CONTR.)
C          ISP = NSPSPZ in block 12, defined for each ICHORI      
C          PSIG(ISP=1), ARGST(ISP=1,) - electron temerature
C          PSIG(2), ARGST (2,) - main ion temperature
C          PSIG(3), ARGST (3,) - electron density
C          PSIG(4), ARGST (4,) - main ion density      
C      
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CGEOM

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: INIT,JJJ
      REAL(DP), INTENT(INOUT) :: PSIG(0:),ARGST(0:,:)
      REAL(DP), INTENT(IN) :: ZDS,DUM1,DUM2
      INTEGER :: NCELC,ITROLD
      LOGICAL :: LARGST
      SAVE

      LARGST = SIZE(ARGST,2) >= NRAD

      IF (INIT.EQ.0) THEN
        PSIG=0.
        IF (LARGST) ARGST=0.

        IF (IITER .NE. ITROLD) THEN
          WRITE (IUNOUT,*) ' sigpla '
           
        ENDIF
        ITROLD=IITER
        RETURN
      ENDIF
 

C  LINE INTEGRAL: PLASMA PARAMETER * CM
C
      ncelc=ncltal(ncell)
      PSIG(1)=PSIG(1)+ZDS*TEIN(NCELC)
      PSIG(2)=PSIG(2)+ZDS*TIIN(1,NCELC)
      PSIG(3)=PSIG(3)+ZDS*DEIN(NCELC)
      PSIG(4)=PSIG(4)+ZDS*DIIN(1,NCELC)

      IF (LARGST) THEN
        ARGST(1,JJJ)=TEIN(NCELC)
        ARGST(2,JJJ)=TIIN(1,NCELC)
        ARGST(3,JJJ)=DEIN(NCELC)
        ARGST(4,JJJ)=DIIN(1,NCELC)
      ENDIF
      
      RETURN
      END
