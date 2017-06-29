c-------------------------------------------------------------------------
      SUBROUTINE GRLN(XX,YY,M)
      USE EIRMOD_PRECISION, ONLY: SP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M
      REAL(SP), DIMENSION(M), INTENT(OUT) :: XX, YY
      XX(1)=0.0
      YY(1)=0.0
      RETURN
      END
