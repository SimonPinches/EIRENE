      module eirmod_cvarusr

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      private

      real*8, allocatable, public :: dVelPrl_dt(:), dVelPerp_dt(:)
      real*8, allocatable, public :: df_dChiPrl(:), dg_dChiPrl(:)
      real*8, allocatable, public :: dg_dChiPerp(:), nue(:)
      real*8, public :: veltotal, old01, old02, old03

      public :: eirene_alloc_cvarusr, eirene_dealloc_cvarusr



      contains

      subroutine eirene_alloc_cvarusr(ical)

      implicit none
      integer, intent(in) :: ical

      if (allocated(dVelPrl_dt)) return

      if (ical == 1) then

         allocate(dVelPrl_dt(1:NPLS))
         allocate(dVelPerp_dt(1:NPLS))
         allocate(df_dChiPrl(1:NPLS))
         allocate(dg_dChiPrl(1:NPLS))
         allocate(dg_dChiPerp(1:NPLS))
         allocate(nue(1:NPLS))
         old01 = 1.0E-10
         old02 = 1.0E-10
         old03 = 1.0E-10

      else if (ical == 2) then

      end if

      call eirene_init_cvarusr(ical)

      return
      end subroutine eirene_alloc_cvarusr



      subroutine eirene_dealloc_cvarusr

      implicit none

         deallocate(dVelPrl_dt)
         deallocate(dVelPerp_dt)
         deallocate(df_dChiPrl)
         deallocate(dg_dChiPrl)
         deallocate(dg_dChiPerp)
         deallocate(nue)

      return
      end subroutine eirene_dealloc_cvarusr


      subroutine eirene_init_cvarusr (ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then

         dVelPrl_dt  = 0.0
         dVelPerp_dt = 0.0
         df_dChiPrl  = 0.0
         dg_dChiPrl  = 0.0
         dg_dChiPerp = 0.0
         nue = 0.0


      else if (ical == 2) then

      end if

      return
      end subroutine eirene_init_cvarusr


      end module eirmod_cvarusr
