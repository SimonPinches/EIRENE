c------------------------------------------------------------------------
      SUBROUTINE GR3EXT(AR,IER,EXT)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: IER
      REAL, INTENT(IN OUT) :: AR(46*128*128), EXT(3,3)
      IER=0
      RETURN
      END SUBROUTINE GR3EXT
