C
C
      SUBROUTINE EIRENE_EXIT_OWN (ICC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ICC
C  LEGAL ENDING OF  GR PLOTTING SOFTWARE      
      CALL EIRENE_PLEND
!pb compiler complains about ICC
!pb no variables allowed, only character or integer constants
!     STOP ICC
      STOP 
      END SUBROUTINE EIRENE_EXIT_OWN
