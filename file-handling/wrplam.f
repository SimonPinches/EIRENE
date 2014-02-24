C  write plasma (background) data, source distribution and atomic data
C  on unit 13.
C
c  at entry RPLAM:
C  read plasma (background) data, source distribution and atomic data
C  from unit 13.
C
C  trcfle:  confirm writing on printout on unit IUNOUT
C  IFLG  :
 
      SUBROUTINE EIRENE_WRPLAM(TRCFLE,IFLG)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE

      IF (NLSHRT13) THEN
        CALL EIRENE_WRPLAM_SHRT (TRCFLE,IFLG)
      ELSE 
        CALL EIRENE_WRPLAM_LONG (TRCFLE,IFLG)
      ENDIF
      RETURN
C
      ENTRY EIRENE_RPLAM(TRCFLE,IFLG)
      IF (NLSHRT13) THEN
        CALL EIRENE_RPLAM_SHRT (TRCFLE,IFLG)
      ELSE 
        CALL EIRENE_RPLAM_LONG (TRCFLE,IFLG)
      ENDIF
      RETURN
      END
