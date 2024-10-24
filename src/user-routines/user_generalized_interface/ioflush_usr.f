      subroutine ioflush_usr
      use eirmod_comprt, only: iunout
      implicit none

      flush(iunout)

      end subroutine ioflush_usr
