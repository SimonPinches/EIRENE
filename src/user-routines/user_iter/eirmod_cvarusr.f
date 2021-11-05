      module eirmod_cvarusr

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      private

      public :: eirene_alloc_cvarusr, eirene_dealloc_cvarusr

      character(15), public, save :: cfort30=repeat(' ',15)
      integer, public, save :: NADMOD=0, NASMOD=0, NORMOD=0

      contains

      subroutine eirene_alloc_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

      else if (ical == 2) then

      end if

      call eirene_init_cvarusr(ical)

      return
      end subroutine eirene_alloc_cvarusr



      subroutine eirene_dealloc_cvarusr

      implicit none

      return
      end subroutine eirene_dealloc_cvarusr


      subroutine eirene_init_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

      else if (ical == 2) then

      end if

      return
      end subroutine eirene_init_cvarusr


      end module eirmod_cvarusr
