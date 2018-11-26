      subroutine usr_defaults

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comprt,only: IUNIN

      implicit none

c     Reset IUNIN to 1
      IUNIN = 1
c     set offset for file units
      IFOFF = 500

      return
      end subroutine usr_defaults
