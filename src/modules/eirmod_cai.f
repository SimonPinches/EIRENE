C> \brief Run time optimisation storage
C>
C> This module stores all variables that are required to allow an 
C> optimisation of the run time based on some simple quantities of a 
C> previous run.
      MODULE EIRMOD_CAI
 
      USE EIRMOD_PRECISION, ONLY: DP
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CAI, EIRENE_DEALLOC_CAI, EIRENE_INIT_CAI
 
C> Ratio between used and recommended no. of particles.
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: RATIO(:)
 
C> Recommended number of test particles for next MC cycle.
      INTEGER, PUBLIC, ALLOCATABLE, SAVE :: NRECOM(:)
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_CAI
 
      USE EIRMOD_PARMMOD, ONLY: IFOFF, NSTRA

      IF (ALLOCATED(RATIO)) RETURN
 
      ALLOCATE (RATIO(NSTRA))
      ALLOCATE (NRECOM(NSTRA))
 
      WRITE (55+IFOFF,'(A,T25,I15)')
     .      ' CAI ',NSTRA*(8+4)
 
      CALL EIRENE_INIT_CAI
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_CAI
 
 
      SUBROUTINE EIRENE_DEALLOC_CAI
 
      IF (.NOT.ALLOCATED(RATIO)) RETURN
 
      DEALLOCATE (RATIO)
      DEALLOCATE (NRECOM)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CAI
 
 
      SUBROUTINE EIRENE_INIT_CAI
 
      RATIO  = 0.D0
      NRECOM = 0
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CAI
 
      END MODULE EIRMOD_CAI
