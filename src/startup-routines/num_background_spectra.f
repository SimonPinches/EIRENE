
      SUBROUTINE EIRENE_NUM_BACKGROUND_SPECTRA(IADTYP)

!     determine number of background spectra

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR, ONLY : ISPZ_BACK, LSPCCLL

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IADTYP(0:4)
      INTEGER :: J, ISPZ

c  number of spectra directly estimated from Monte Carlo trajectories

      NADSPC_S = 0   !  surface-based
cdr   NADSPC_P = 0   !  surface-based directional missing? Already programmed in OUTSPEC?
cdr   NADSPC_SP = 0  !  surface-based, both
      NADSPC_C = 0   !  cell-based
      NADSPC_D = 0   !  cell-based, directional
      NADSPC_CD = 0  !  cell-based, both

!   Determine these numbers of spectra (from block 10F), by type,
cdr as well as the number NBACK_SPEC of particular directional cell based
cdr spectra needed for line of sight integration
cdr as specified by an option in input block 12.

      NBACK_SPEC = 0
      DO J = 1, NADSPC
!  directional spectrum in geometrical cell
cdr I do not know what that option is. Unfinished or redundant?
cdr Partially used in OUTSPEC, but not for calling UPDATE_SPECTRUM.
        IF (ESTIML(J)%ISRFCLL == 2) THEN
          ISPZ=IADTYP(ESTIML(J)%IPRTYP) + ESTIML(J)%IPRSP
          NBACK_SPEC = NBACK_SPEC + COUNT(ISPZ_BACK(ISPZ,:)>0)
          LSPCCLL(ESTIML(J)%ISPCSRF) = .TRUE.
        END IF

        IF (ESTIML(J)%ISRFCLL == 0) THEN
C  COUNT SURFACE SPECTRA
          NADSPC_S=NADSPC_S+1
cdr  later: idirec > 0 is used to distinguish total from directional spectra.
cdr  here we should add directional surface spectra. In OUTSPEC already done.

        ELSEIF (ESTIML(J)%ISRFCLL == 1) THEN
cdr  later: isrfcll seems to have a different meaning.
cdr         =1: score in scoring cell (of coarse grid)
cdr         =2  score in geometry cell (of fine grid)
C  COUNT CELL-BASED SPECTRA
          NADSPC_C=NADSPC_C+1
        ELSEIF (ESTIML(J)%ISRFCLL == 2) THEN
C  COUNT DIRECTIONAL CELL-BASED SPECTRA
cdr  in other parts of the code isrfcll seems to
cdr  switch between (fine) geometry grid cells and (coarse) scoring grid cells
cdr  whereas IDIREC is used to switch between cell averages and directional spectra
          NADSPC_D=NADSPC_D+1
        ENDIF
      END DO
      NADSPC_CD=NADSPC_C+NADSPC_D

      RETURN
      END SUBROUTINE EIRENE_NUM_BACKGROUND_SPECTRA
