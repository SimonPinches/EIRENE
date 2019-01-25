      subroutine eirene_defaults_usr

      use eirmod_precision
      use eirmod_parmmod

      implicit none

c     set file units that allow direct control
      IUNMEM = 85
      IUNRAPS = 80
      IUNRAPSVEC = 70

      return
      end subroutine eirene_defaults_usr
