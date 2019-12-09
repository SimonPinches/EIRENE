c------------------------------------------------------------------------
      SUBROUTINE GR3NT1(AR,IER,I1,X,Y,Z,I2,I3,I4,I5,I6,I7)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IER
      INTEGER, INTENT(IN) :: I1, I2, I3, I4, I5, I6, I7
      REAL, INTENT(IN OUT) :: AR(46*I1*1)
      REAL, INTENT(IN) :: X(I1,1), Y(I1,1), Z(I1,1)
      IER=0
      RETURN
      END
