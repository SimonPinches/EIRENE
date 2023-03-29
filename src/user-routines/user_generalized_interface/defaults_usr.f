      subroutine eirene_defaults_usr

      use eirmod_parmmod
      use eirmod_cpes, only: STRATEGY_DEFAULT, NPRLL_DEFAULT,
     .                       STRATEGY_EMBARRASS

      implicit none

      loutapp = .false.
      lif3cop_from_loop = .false.

      STRATEGY_DEFAULT = STRATEGY_EMBARRASS
      NPRLL_DEFAULT = 0

      return
      end subroutine eirene_defaults_usr
