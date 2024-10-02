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
      SUBROUTINE EIRENE_WRPLAM(TRCFLE,CALLEDFROM)
cdr called from subr. 'calledfrom'
C  trcfle:  confirm writing on printout unit IUNOUT

cdr  NLSHRT13 : VIA COMMON CLOGAU. MEANING: write "long" or "short" version of fort.13
C               DEFAULT: FALSE,
cdr  but, e.g.:
cdr  NLSHRT13 : SET TRUE IN INFCOP, COUPLE_SOLPS_ITER. REDUCED SIZE FORT.13.
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU, ONLY: NLSHRT13

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: TRCFLE
      CHARACTER(*), INTENT(IN) :: CALLEDFROM

      IF (NLSHRT13) THEN
        CALL EIRENE_WRPLAM_SHRT (TRCFLE,CALLEDFROM)
      ELSE
        CALL EIRENE_WRPLAM_LONG (TRCFLE,CALLEDFROM)
      ENDIF
      RETURN
      END SUBROUTINE EIRENE_WRPLAM
C
c.............................................

      SUBROUTINE EIRENE_RPLAM(TRCFLE,IFLG,CALLEDFROM)
cdr called from subr. 'calledfrom'
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
      USE EIRMOD_CLOGAU, ONLY: NLSHRT13
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEI
      USE EIRMOD_CINIT

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IFLG
      LOGICAL, INTENT(IN) :: TRCFLE
      CHARACTER(*), INTENT(IN) :: CALLEDFROM
      INTEGER :: JPLS
c.............................................

      IF (NLSHRT13) THEN  !dr  similar to RPLAM_LONG(IFLG=10), but for
cdr                            field species  NLFA+1:NPLS only.
        CALL EIRENE_RPLAM_SHRT(TRCFLE,CALLEDFROM)

cdr  only species NPLS_FIX+1:NPLS should be affected from RPLAM_SHRT.
cdr  No magnetic field etc..
        CALL EIRENE_ALLOC_BCKGRND
        TEINTF(1:NRAD) = TEIN(1:NRAD)
        IF (NLMLTI) THEN
          DO JPLS = 1, NPLSI
            TIINTF(MPLSTI(JPLS),1:NRAD) = TIIN(JPLS,1:NRAD)
          ENDDO
        ELSE
          TIINTF(1,1:NRAD) = TIIN(1,1:NRAD)
        ENDIF
        DIINTF(1:NPLSI,1:NRAD) = DIIN(1:NPLSI,1:NRAD)
        IF (NLMLV) THEN
          DO JPLS = 1, NPLSI
            VXINTF(MPLSV(JPLS),1:NRAD) = VXIN(JPLS,1:NRAD)
            VYINTF(MPLSV(JPLS),1:NRAD) = VYIN(JPLS,1:NRAD)
            VZINTF(MPLSV(JPLS),1:NRAD) = VZIN(JPLS,1:NRAD)
          ENDDO
        ELSE
          VXINTF(1,1:NRAD) = VXIN(1,1:NRAD)
          VYINTF(1,1:NRAD) = VYIN(1,1:NRAD)
          VZINTF(1,1:NRAD) = VZIN(1,1:NRAD)
        ENDIF
      ELSE
cdr  set pointers teintf,.... only in case iflg=10  (same as nlshrt13?)
        CALL EIRENE_RPLAM_LONG (TRCFLE,IFLG,CALLEDFROM)
      ENDIF
      RETURN

      END SUBROUTINE EIRENE_RPLAM
