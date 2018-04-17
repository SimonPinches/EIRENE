*DK INFCOP
      SUBROUTINE EIRENE_INFCOP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1, I2, I3
      ENTRY EIRENE_IF0COP
      ENTRY EIRENE_IF1COP
      ENTRY EIRENE_IF2COP(I1)
      ENTRY EIRENE_IF3COP(I1,I2,I3)
      ENTRY EIRENE_IF4COP
      END

C> \brief Any property requirering hand-over in parallel part.
C>
C> This interfacing routine is called in the parallel part of EIRENE
C> after the broadcase of any other quantity.
      SUBROUTINE EIRENE_IFPARCOP
      RETURN
      END SUBROUTINE EIRENE_IFPARCOP
