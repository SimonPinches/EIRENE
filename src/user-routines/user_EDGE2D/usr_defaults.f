      subroutine eirene_defaults_usr

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comprt,only: IUNIN

      implicit none

c     Reset IUNIN to 1
      IUNIN = 1
c     set offset for file units
      IFOFF = 500
c     set file units that allow direct control
      IUNMEM = 555
      IUNRAPS = 565
      IUNRAPSVEC = 560

      return
      end subroutine eirene_defaults_usr
