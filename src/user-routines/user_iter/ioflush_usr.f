      subroutine ioflush_usr
      implicit none

#ifdef WINDOWS
      interface
      subroutine ioflush
!DIR$attributes c, alias: 'ioflush_' :: ioflush
      end subroutine
      end interface
#endif

cpg ioflush missing on julich repo
cpg        call ioflush

      end subroutine ioflush_usr
