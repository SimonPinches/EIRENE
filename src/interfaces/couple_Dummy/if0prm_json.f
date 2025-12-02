cdr called from find_param.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C
#ifndef NO_JSON
      SUBROUTINE EIRENE_IF0PRM_JSON(json,p)

      use eirmod_parmmod, only: ncpv, nain, nptrgt
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      IMPLICIT NONE

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      NCPV=0
      NAIN=0
      NPTRGT=1

      RETURN
      END SUBROUTINE EIRENE_IF0PRM_JSON
#else
      SUBROUTINE EIRENE_IF0PRM_JSON
      RETURN
      END SUBROUTINE EIRENE_IF0PRM_JSON
#endif
