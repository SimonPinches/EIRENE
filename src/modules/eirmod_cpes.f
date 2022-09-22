cdr Nov. 17  commenting started
cdr July 18  remove nsteff, redundant

      MODULE EIRMOD_CPES

      USE EIRMOD_PARMMOD, ONLY: IUNMEM, NSTRA
      USE EIRMOD_MPI, ONLY: MPI_COMM_NULL, MPI_SET_OWN_IO_UNIT
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CPES, EIRENE_DEALLOC_CPES, EIRENE_INIT_CPES
      public :: I_am_leader
      public :: create_all_communicators
      public :: get_leader_comm, get_stratum_comm
      public :: need_calstr, calc_stratum

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
cdr  npesta(istra): no. of master processor for ISTRA
cdr  npestr(istra): total number of processors working on ISTRA
     I         NPESTR(:), NPESTA(:)

      INTEGER, PUBLIC, SAVE :: NPRS, MY_PE

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     I         LEXIT(:)

      LOGICAL, PUBLIC, SAVE :: NLIDENT

CVKMPI CORRESPONDENCE TABLE "STRATA VERSUS PROCESSOR"
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: PROCFORSTRA(:,:)

      integer, dimension(:), allocatable, public :: stratum_comm 
      integer, public :: leader_comm = MPI_COMM_NULL

      !> Whether we have new member of stratum leaders
      logical, public :: new_leader = .true. 

      CONTAINS

      !> creates communicators for PEs working on same stratum and a 
      !> communicator for all stratum leaders
      subroutine create_all_communicators()
        use eirmod_comsou, only: nlsron
        use eirmod_mpi, only: MPI_COMM_NULL
        use eirmod_parmmod, only: nstra
        integer :: k
        call free_communicators(nstra)
        ! Initialize communicators 
        do k=1, nstra
          if (nlsron(k)) then
            stratum_comm(k) = create_stratum_comm(k)
          else
            stratum_comm(k) = MPI_COMM_NULL
          endif
        end do
        leader_comm = create_leader_comm()
      end subroutine

      !> returns true if the calling PE should do any work on stratum_idx
      logical function calc_stratum(stratum_idx)
        integer, intent(in) :: stratum_idx
        calc_stratum = procforstra(stratum_idx, my_pe)
      end function

      !> creates the communicator for all stratum leaders
      function create_leader_comm() result(comm)
        integer :: comm
        if (my_pe==0) then
          write(iunout,*) 'Creating communicator for stratum leaders'
        end if
        comm = create_communicator(I_am_leader())
      end function

      !> creates the communicators for PEs working on same stratum
      function create_stratum_comm(stratum_idx) result(comm)
        integer :: comm
        integer, intent(in) :: stratum_idx
        if (my_pe==0) then
          write(iunout,*) 'Creating communicator within stratum ',
     &     stratum_idx
        end if
        comm = create_communicator(calc_stratum(stratum_idx), 
     &           npesta(stratum_idx))
      end function      

      !> creates a communicator for the PEs which call it with .true. argument
      !> It is a collective call: every PE in MPI_COMM_WORLD should call it together.
      function create_communicator(include_me, zero_rank) result(comm)
        use eirmod_mpi
        !> whether to include the calling PE
        logical, intent(in) :: include_me 
        !> the rank of the PE which should become rank 0 in the new communicator
        integer, optional :: zero_rank 
        integer :: comm 
        integer :: n, k, ierr
        integer, dimension(nprs) :: ranks_tmp
        logical, dimension(nprs) :: include_proc
        integer :: group_world
        integer :: group !, gsize, grank
        logical :: found
        character*13 hlp_frm

        call MPI_allGather(include_me, 1, MPI_LOGICAL, 
     &          include_proc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
        ! now collect all the ranks where include_proc is true
        n = 0
        do k=0, nprs-1
          if (include_proc(k+1)) then
            n = n + 1
            ranks_tmp(n) = k
          end if
        end do

        if (present(zero_rank)) then
          ! put zero_rank into ranks_tmp(1)
          ! but first free the place ranks_tmp(1)
          found = .false.
          do k=1,n
            if (ranks_tmp(k) == zero_rank) then
              ranks_tmp(k) = ranks_tmp(1)
              found = .true.
            endif
          end do
          if (found) then
            ranks_tmp(1) = zero_rank
          else
            write(iunout,*) "Error, rank ", zero_rank,
     &         " is not a member of the new communicator"
          endif
        end if
        if (my_pe == 0) then
          write(hlp_frm,'(a,i4,a)') '(1x,a,',n,'i5)'
          write(iunout,hlp_frm)
     &     'Creating group from ranks ', ranks_tmp(1:n)
        end if

        call mpi_comm_group(MPI_COMM_WORLD, group_world, ierr)
        call mpi_group_incl(group_world, n, ranks_tmp, group, ierr)
        call mpi_comm_create(MPI_COMM_WORLD, group, comm, ierr)
        !if (include_me) then
        !  call mpi_group_size(group, gsize, ierr)
        !  call mpi_group_rank(group, grank, ierr)
        !  write(iunout,*) 'My rank in the new communicator ', grank, ' (', gsize,')'
        !endif
        call mpi_group_free(group_world, ierr)
        call mpi_group_free(group, ierr)
      end function

      !> frees all communicators
      subroutine free_communicators(nstra)
        use eirmod_mpi
        integer, intent(in) ::nstra
        integer k, ierr
        do k = 1, nstra
          if(stratum_comm(k) .ne. MPI_COMM_NULL) then
            call mpi_comm_free(stratum_comm(k), ierr)
          end if
        end do
        if (leader_comm .ne. MPI_COMM_NULL) then
          call mpi_comm_free(leader_comm, ierr)
        endif
      end subroutine

      !> returns the communicator of the stratum leader
      function get_leader_comm() result(comm)
        integer :: comm
        comm = leader_comm
      end function

      !> returns the communicator of the PEs working on stratum `stratum_idx`
      function get_stratum_comm(stratum_idx) result(comm)
        integer, intent(in) :: stratum_idx
        integer :: comm
        comm = stratum_comm(stratum_idx)
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

      SUBROUTINE EIRENE_ALLOC_CPES(ICAL)

      INTEGER ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(NPESTR)) RETURN

        ALLOCATE (NPESTR(NSTRA))
        ALLOCATE (NPESTA(NSTRA))

        ALLOCATE(PROCFORSTRA(NSTRA,0:NPRS-1))

        ALLOCATE(STRATUM_COMM(NSTRA))
        STRATUM_COMM = MPI_COMM_NULL

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' CPES(1) ',2*NSTRA*4 + NSTRA*NPRS*4

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(LEXIT)) RETURN
        ALLOCATE (LEXIT(0:NPRS-1))

        WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CPES(2) ', NPRS*4

      END IF

      CALL EIRENE_INIT_CPES(ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CPES


      SUBROUTINE EIRENE_DEALLOC_CPES

      IF (.NOT.ALLOCATED(NPESTR)) RETURN

      DEALLOCATE (NPESTR)
      DEALLOCATE (NPESTA)

      DEALLOCATE (LEXIT)

      DEALLOCATE(PROCFORSTRA)

      CALL FREE_COMMUNICATORS(NSTRA)
      DEALLOCATE(STRATUM_COMM)
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CPES


      SUBROUTINE EIRENE_INIT_CPES(ICAL)
      INTEGER ICAL

      IF (ICAL == 1) THEN
      
        NPESTR = 0
        NPESTA = 0

c  correspondence table: Strata vs. PEs
        PROCFORSTRA=.TRUE.      !  Trivial parallelisation: All PEs work on all strata

      ELSE IF (ICAL == 2) THEN

        LEXIT = .FALSE.

      END IF

      CALL MPI_SET_OWN_IO_UNIT(IUNOUT)
      RETURN
      END SUBROUTINE EIRENE_INIT_CPES

      END MODULE EIRMOD_CPES
