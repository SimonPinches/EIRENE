      module eirmod_cvarusr

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      private

      real*8, allocatable, public :: dVelPrl_dt(:), dVelPerp_dt(:)

      public :: eirene_alloc_cvarusr, eirene_dealloc_cvarusr



      contains

      subroutine eirene_alloc_cvarusr(ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

         allocate(dVelPrl_dt(1:NPLS))
         allocate(dVelPerp_dt(1:NPLS))

      else if (ical == 2) then

      end if

      call eirene_init_cvarusr(ical)

      return
      end subroutine eirene_alloc_cvarusr



      subroutine eirene_dealloc_cvarusr

      implicit none

         deallocate(dVelPrl_dt)
         deallocate(dVelPerp_dt)

      return
      end subroutine eirene_dealloc_cvarusr


      subroutine eirene_init_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

         dVelPrl_dt  = 0.0
         dVelPerp_dt = 0.0

      else if (ical == 2) then

      end if

      return
      end subroutine eirene_init_cvarusr


      end module eirmod_cvarusr
