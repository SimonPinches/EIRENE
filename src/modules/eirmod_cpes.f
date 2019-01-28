cdr Nov. 17  commenting started
cdr July 18  remove nsteff, redundant

      MODULE EIRMOD_CPES

      USE EIRMOD_PARMMOD, ONLY: IUNMEM, NSTRA

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CPES, EIRENE_DEALLOC_CPES, EIRENE_INIT_CPES
      public :: I_am_leader
      public :: need_calstr, calc_stratum

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
cdr  npesta(istra): master processor for ISTRA
cdr  npestr(istra): total no. of processor working on ISTRA
     I         NPESTR(:), NPESTA(:)

      INTEGER, PUBLIC, SAVE :: NPRS, MY_PE

      LOGICAL, PUBLIC, SAVE :: NLIDENT

CVKMPI CORRESPONDENCE TABLE "STRATA VERSUS PROCESSOR"
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: PROCFORSTRA(:,:)


      CONTAINS

      !> returns true if the calling PE should do any work on stratum_idx
      logical function calc_stratum(stratum_idx)
        integer, intent(in) :: stratum_idx
        calc_stratum = procforstra(stratum_idx, my_pe)
      end function

      !> Checks whether the calling PE is the leader of stratum with stratum_idx,
      !> or any stratum if stratum_idx is not present
      logical function I_am_leader(stratum_idx)
        use eirmod_comsou, only: nlsron
        integer, optional, intent(in) :: stratum_idx
        if (present(stratum_idx)) then
          I_am_leader = npesta(stratum_idx) == my_pe
        else
          I_am_leader = any(npesta==my_pe .and. nlsron)
        end if
      end function

      !> returns true if the calling PE should call eirene_calstr
      logical function need_calstr(stratum_idx)
        integer, intent(in) :: stratum_idx
        need_calstr = npestr(stratum_idx) > 1
      end function

      SUBROUTINE EIRENE_ALLOC_CPES

      IF (ALLOCATED(NPESTR)) RETURN

      ALLOCATE (NPESTR(NSTRA))
      ALLOCATE (NPESTA(NSTRA))

      ALLOCATE(PROCFORSTRA(NSTRA,0:NPRS-1))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CPES ',2*NSTRA*4 + NSTRA*NPRS*4

      CALL EIRENE_INIT_CPES

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CPES


      SUBROUTINE EIRENE_DEALLOC_CPES

      IF (.NOT.ALLOCATED(NPESTR)) RETURN

      DEALLOCATE (NPESTR)
      DEALLOCATE (NPESTA)

      DEALLOCATE(PROCFORSTRA)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CPES


      SUBROUTINE EIRENE_INIT_CPES

      NPESTR = 0
      NPESTA = 0

      PROCFORSTRA=.TRUE. !VK

      RETURN
      END SUBROUTINE EIRENE_INIT_CPES

      END MODULE EIRMOD_CPES
