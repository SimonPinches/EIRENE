      subroutine eirene_wr0_json(json,me)

      use eirmod_precision
      use json_module           !IGNORE

      implicit none

      class(json_core),intent(inout) :: json
      type(json_value),pointer, intent(inout) :: me

      end subroutine eirene_wr0_json
