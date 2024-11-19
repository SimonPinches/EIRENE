c-------------------------------------------------------------------------
      SUBROUTINE GRFILL(N,XX,YY,ISTYLE,ITYPE)
      USE EIRMOD_PRECISION, ONLY: SP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, ISTYLE, ITYPE
      REAL(SP), DIMENSION(N), INTENT(IN) :: XX, YY
      RETURN
      END SUBROUTINE GRFILL
