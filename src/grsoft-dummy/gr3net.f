c------------------------------------------------------------------------
      SUBROUTINE GR3NET(AR,IER,I1,XYZ,I2,I3,I4,I5,I6,I7)
      IMPLICIT NONE
      REAL, INTENT(OUT) :: AR(*)
      REAL, INTENT(IN) :: XYZ(3,*,*)
      INTEGER, INTENT(IN) :: IER, I1, I2, I3, I4, I5, I6, I7
      AR = 0
      RETURN
      END
