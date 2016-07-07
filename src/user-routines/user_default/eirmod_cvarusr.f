      module eirmod_cvarusr

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IION
      USE EIRMOD_COMXS, ONLY: NIELI

      IMPLICIT NONE

      private

      real(dp), allocatable, public :: dVelPrl_dt(:), dVelPerp_dt(:)
      real(dp), allocatable, public :: df_dChiPrl(:), dg_dChiPrl(:)
      real(dp), allocatable, public :: dg_dChiPerp(:), nue(:)
      real(dp), public  :: veltotal, rCPrlOld, rCPerpOld
C     real(dp), public  :: old03, old04, old05, old06, old07
      real(dp), public  :: alphaPrl, alphaPerp, iprepare
      integer, public :: npanuSave

      public :: eirene_alloc_cvarusr, eirene_dealloc_cvarusr



      contains

      subroutine eirene_alloc_cvarusr(ical)

      implicit none
      integer, intent(in) :: ical

      if (ical == 1) then
        if (.not.allocated(dVelPrl_dt)) then
          allocate(dVelPrl_dt(1:NIELI(IION)))
          allocate(dVelPerp_dt(1:NIELI(IION)))
          allocate(df_dChiPrl(1:NIELI(IION)))
          allocate(dg_dChiPrl(1:NIELI(IION)))
          allocate(dg_dChiPerp(1:NIELI(IION)))
          allocate(nue(1:NIELI(IION)))
          rCPrlOld  = 1.0E-10
          rCPerpOld = 1.0E-10
C         old03 = 1.0E-10
C         old04 = 1.0E-10
C         old05 = 1.0E-10
C         old06 = 1.0E-10
C         old07 = 1.0E-10
          npanuSave = 0
          iprepare = 0.0
        end if
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
         nue = 1.D-30

      else if (ical == 2) then

      end if

      return
      end subroutine eirene_init_cvarusr


      end module eirmod_cvarusr
