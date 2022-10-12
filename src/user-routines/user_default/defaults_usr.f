      subroutine eirene_defaults_usr

      use eirmod_parmmod
      use eirmod_cpes, only: nprs

      implicit none

      loutapp = .false.
      lif3cop_from_loop = .false.
      if (nprs == 1) lif3cop_from_loop = .true.

      return
      end subroutine eirene_defaults_usr
