cdr called from find_param.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C
      SUBROUTINE EIRENE_IF0PRM_JSON(json,p)

      use json_module
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      IMPLICIT NONE

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      NCPV=0
      NAIN=0
      NPTRGT=1

      RETURN
      END
