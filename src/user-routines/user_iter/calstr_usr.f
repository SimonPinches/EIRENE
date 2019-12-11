
      SUBROUTINE EIRENE_CALSTR_USR(MY_PE, ICGRP)
c called from broadcast/calstr
c can be used from user or case-specific routines (...usr.f, ...cop.f)
c to collect user/case  specific data from all
c PEs sharing calculations for a particular stratum
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MY_PE, ICGRP
      RETURN
      END
