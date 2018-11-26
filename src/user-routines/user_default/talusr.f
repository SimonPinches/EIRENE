c
c
      subroutine EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(inout) :: VECTOR(*), TALTOT, TALAV

      character(len=*) :: txttl,txtsp,txtun
      TXTTL=' '
      TXTSP=' '
      TXTUN=' '


      ILAST=0
      return
      end
