cdr
c  at subroutine WRPLAM:
C  write plasma (background) data, primary source distribution and atomic data
C  onto unit fort.13.
C
cdr
c  at subroutine RPLAM:
C  read plasma (background) data, primary source distribution and atomic data
C  from unit 13. Allocate target: plasma_bckgrnd

c  if NLSHRT13=true : wrplam_short and rplam_short are called from
c                     WRPLAM, RPLAM, resp.
c                     only virtual background (BGK-species) NPLS_FIX+1:NPLS data
c                     are stored on on fort.13.
c                     No primary source data
c                     Only special subsets of atomic rates and related info.
c                     are on fort.13:
c  if NLSHRT13=false: wrplam_long and rplam_long are called from
c                     WRPLAM, RPLAM, resp.
c
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

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL, INTENT(IN) :: TRCFLE

      IF (NLSHRT13) THEN
        CALL EIRENE_WRPLAM_SHRT
      ELSE
        CALL EIRENE_WRPLAM_LONG (TRCFLE,IFLG)
      ENDIF
      END SUBROUTINE EIRENE_WRPLAM
C
c.............................................

      SUBROUTINE EIRENE_RPLAM(TRCFLE,IFLG,IRET)
C  TRCFLE:    confirm writing on printout on unit IUNOUT
C  NLSHRT13: =T: reading only for virt. background species NPLS_FIX+1,NPLS:
C                call rplam_shrt, allocate plasma_bckgrnd,
C                and set pointers TEINTF,... for virt. species NPLS_FIX+1,NPLS
C            =F: call rplam_long, full COMUSR, A&M DATA, COMSOU
cdr
C                IFLG:  only for RPLAM_LONG version
C                fort.13 was written in the "long" format.
C                IFLG = 0  : read full COMUSR, A&M DATA, COMSOU
c                IFLG > 0  : same, but do not read COMSOU as well.
c                IFLG =10  : something special: do not read A&M data, nor primary source data,
cdr                          and read only the input field particle tally part PLSTLS from COMUSR
cdr                          but no other parameters.
cdr                          Allocate plasma_bckgrnd, and set pointers Teintf....
cdr                          This option is similar to the NLSHRT13 mode, but it handles
cdr                          all ipls: 1,npls, not just ipls=nfla+1,npls

      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEI

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      INTEGER, INTENT(INOUT) :: IRET
      LOGICAL, INTENT(IN) :: TRCFLE
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

      END SUBROUTINE EIRENE_RPLAM
