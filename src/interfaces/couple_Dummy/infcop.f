*DK INFCOP
      MODULE EIRMOD_INFCOP
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE EIRENE_INFCOP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I1, I2, I3
      ENTRY EIRENE_IF0COP
      ENTRY EIRENE_IF1COP
      ENTRY EIRENE_IF2COP(I1)
      ENTRY EIRENE_IF3COP(I1,I2,I3)
      ENTRY EIRENE_IF4COP
      END SUBROUTINE

C> \brief Any property requiring hand-over in parallel part.
C>
C> This interfacing routine is called in the parallel part of EIRENE
C> after the broadcast of any other quantity and before MCARLO.
      SUBROUTINE EIRENE_INFCOP_PRE_MCARLO
      RETURN
      END SUBROUTINE EIRENE_INFCOP_PRE_MCARLO

C> \brief Return data at the end of an EIRENE stratum
C>
C> At this entry data are transferred back from EIRENE to the
C> interfacing module (and from there, after possibly further
C> processing, to the external code). This entry is called from the
C> "strata-loop" in subr. MCARLO, after all trajectories for a
C> particular stratum ISTRA have been sampled and after all volume and
C> surface tallies have been scaled and processed to their final form.
C>
C> The call to IF3COP is controlled by the flag NMODE (input block 1).
C> At call data are expected for stratum ISTRA. There they may be
C> further prepared (e.g. normalized, or scaled to other units) for
C> transfer to the external code
      SUBROUTINE EIRENE_INFCOP_POST_STRATUM(ISTRA)
      RETURN
      END SUBROUTINE EIRENE_INFCOP_POST_STRATUM

      END MODULE
