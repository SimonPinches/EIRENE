
      SUBROUTINE EIRENE_CALSTR_USR(MY_PE, ICGRP)
      ! This subroutine is not called if you use STRATEGY_BALANCED parallelizaton.
      ! In that case only eirene_calstr_buffered is called.
      ! Therefore, if you add any reduction operations here, then please add the same to
      ! eirmod_calstr_buffered.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MY_PE, ICGRP
      RETURN
      END SUBROUTINE EIRENE_CALSTR_USR
