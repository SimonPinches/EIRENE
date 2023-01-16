c------------------------------------------------------------------------
      SUBROUTINE GR3NET(AR,IER,I1,XYZ,I2,I3,I4,I5,I6,I7)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IER, I1, I2, I3, I4, I5, I6, I7
      REAL, INTENT(IN OUT) :: AR(46*I1*I1)
      REAL, INTENT(IN) :: XYZ(3,I1,I1)
      RETURN
      END SUBROUTINE GR3NET
