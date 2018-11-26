C> \brief Initialisation of additional particle property at birth.
C>
C> User routine that allows to set some additional particle properties
C> required by the plasma code.
      SUBROUTINE EIRENE_CNSUSR( NPANU )

      IMPLICIT NONE
C> Particle number
      INTEGER, INTENT(IN) :: NPANU

      RETURN
      END SUBROUTINE EIRENE_CNSUSR
