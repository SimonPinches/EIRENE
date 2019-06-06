cdr
c  at entry WRPLAM:
C  write plasma (background) data, source distribution and atomic data
C  onto unit fort.13.
C
cdr
c  at entry RPLAM:
C  read plasma (background) data, source distribution and atomic data
C  from unit 13.
C
C  trcfle:  confirm writing on printout on unit IUNOUT

C  iflg:  only for .._LONG version, and there only for RPLAM
C  IFLG = 0  :  do NOT read COMSOU in call RPLAM_LONG
c  IFLG > 0  :

cdr  NLSHRT13  :  VIA COMMON CLOGAU. MEANING: write "long vs. short" version of fort.13
C                 DEFAULT: FALSE,
cdr  EXCEPT:
cdr  NLSHRT13  :  SET TRUE IN INFCOP, COUPLE_SOLPS_ITER. REDUCED SIZE FORT 13.

      SUBROUTINE EIRENE_WRPLAM(TRCFLE,IFLG)
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEI

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      INTEGER, INTENT(INOUT) :: IRET
      LOGICAL, INTENT(IN) :: TRCFLE

      IF (NLSHRT13) THEN
        CALL EIRENE_WRPLAM_SHRT (TRCFLE)
      ELSE
        CALL EIRENE_WRPLAM_LONG (TRCFLE,IFLG)
      ENDIF
      RETURN
C
c.............................................

      ENTRY EIRENE_RPLAM(TRCFLE,IFLG,IRET)
c.............................................

      IRET = 0
      IF (NLSHRT13) THEN
        CALL EIRENE_RPLAM_SHRT (TRCFLE)
        CALL EIRENE_ALLOC_BCKGRND
        TEINTF = TEIN
        TIINTF = TIIN
        DIINTF = DIIN
        VXINTF = VXIN
        VYINTF = VYIN
        VZINTF = VZIN
      ELSE
        CALL EIRENE_RPLAM_LONG (TRCFLE,IFLG,IRET)
      ENDIF
      RETURN

      END SUBROUTINE EIRENE_WRPLAM
