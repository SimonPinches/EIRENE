C
C
      SUBROUTINE EIRENE_FTCRE (F,C)
C  write a real number F in format E10.3 onto character string C(10)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: F
      CHARACTER(10), INTENT(OUT) :: C
      WRITE(C,'(1P,E10.3)') F
      RETURN
      END
