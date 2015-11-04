      module eirmod_cvarusr

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      private

      public :: eirene_alloc_cvarusr, eirene_dealloc_cvarusr

      real(DP),allocatable,public,save :: dVelPrl_dt(:), dVelPerp_dt(:)


      contains

      subroutine eirene_alloc_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      IF (allocated(dVelPrl_dt)) RETURN

      allocate(dVelPrl_dt(1:NPLS))
      allocate(dVelPerp_dt(1:NPLS))

      if (ical == 1) then

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

      dVelPrl_dt  = 0._DP
      dVelPerp_dt = 0._DP

      if (ical == 1) then

      else if (ical == 2) then

      end if

      return
      end subroutine eirene_init_cvarusr


      end module eirmod_cvarusr
