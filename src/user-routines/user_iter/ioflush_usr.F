      subroutine ioflush_usr
      use eirmod_comprt, only: iunout
      implicit none

#ifdef WINDOWS
      interface
      subroutine ioflush
!DIR$attributes c, alias: 'ioflush_' :: ioflush
      end subroutine
      end interface
#endif

      flush(iunout)

      end subroutine ioflush_usr
