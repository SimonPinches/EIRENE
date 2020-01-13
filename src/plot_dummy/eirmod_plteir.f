      MODULE EIRMOD_PLTEIR
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_PLTEIR, EIRENE_PLTEIR_REINIT

      CONTAINS

c------------------------------------------------------------------------
      subroutine EIRENE_plteir(i)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      return
      end subroutine EIRENE_plteir

      SUBROUTINE eirene_plteir_reinit
      IMPLICIT NONE
      return
      end SUBROUTINE eirene_plteir_reinit

      end MODULE EIRMOD_PLTEIR
