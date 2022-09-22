      subroutine ioflush_usr
      use eirmod_comprt, only: iunout
      implicit none

      call flush(iunout)

      end subroutine ioflush_usr
