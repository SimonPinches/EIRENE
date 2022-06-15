      subroutine eirene_defaults_usr

      use eirmod_parmmod

      implicit none

c     set file units that allow direct control
      IUNMEM = 85
      IUNRAPS = 80
      IUNRAPSVEC = 70
      LOUTAPP = .FALSE.

      return
      end subroutine eirene_defaults_usr
