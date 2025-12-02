! New balanced strategy optimizing the workload distribution over MPI processes.
! It has been developed based on the original balanced strategy for
! more stable and bettter parallelization performance.
module eirmod_improved_strategy

  use eirmod_precision
  use eirmod_comprt, only: iunout
  use eirmod_comsou, only: npts, nstrai
  use eirmod_parmmod, only: nstra
  implicit none

  private

  ! Init_improved_strategy should be called first, afterwards in every
  ! iteration we optimize
  public init_improved_strategy
  public opt_improved_strategy
  public allocate_improved_strategy
  public deallocate_improved_strategy

  ! The time-measurement subroutines called from mcarlo subroutine.
  ! The workload distribution is optimized based on these data.
  public time_calstr_improved
  public time_postproc_improved
  public time_particles_improved

  ! Distributed throughput and overhead time used only for
  ! printing work distribution table in eirmod_cpes.F.
  ! (allocated only on PE 0)
  real(kind=dp), allocatable, dimension(:), public :: throughput2, t_overhead2

  ! Arrays to store timing information, all local to the current PE:

  ! The n_particles_loc(k) stores how many particles have been
  ! processed on stratum(k) by the local PE.
  ! It is a cumulative quantity, in every iteration we increase
  ! n_particles_loc(k) by the number of particles that were processed.
  integer, allocatable, dimension(:) :: n_particles_loc

  ! t_particles_loc(k) stores the CPU time (in seconds) that was used to
  ! process the particles n_particles_loc(k)
  real(kind=dp), allocatable, dimension(:) :: t_particles_loc, t_particles_prev

  ! How many times calstr was called for each stratum
  integer, allocatable, dimension(:) :: n_calstr_loc

  ! If calstr is called for stratum k, then t_calstr(k) is the execution time
  ! of calstr. This is already a weighted average over iterations of eirene
  real(kind=dp), allocatable, dimension(:) :: t_calstr_loc, t_calstr_prev

  ! The number of all the postprocessing steps so far is
  integer :: n_postproc_loc

  ! t_postproc(k) is the postprocessing time (after calstr) for stratum k
  real(kind=dp), allocatable, dimension(:) :: t_postproc_loc, t_postproc_prev

  ! We cannot acces these from CPES (circular reference) so we store them here
  ! set in init_improved_strategy()
  integer :: my_pe, nprs

  ! Whether to print debug information
#ifdef TRACE
  logical, parameter :: print_debug_info = .true.
#else
  logical, parameter :: print_debug_info = .false.
#endif

  contains

  subroutine allocate_improved_strategy
    allocate(n_particles_loc(nstra))
    allocate(t_particles_loc(nstra))
    allocate(n_calstr_loc(nstra))
    allocate(t_calstr_loc(nstra))
    allocate(t_postproc_loc(nstra))
    call reset_counters()
    if (my_pe == 0) then
      allocate(t_calstr_prev(nstra))
      allocate(t_particles_prev(nstra))
      allocate(t_postproc_prev(nstra))
    endif
    return
  end subroutine allocate_improved_strategy

  subroutine deallocate_improved_strategy
    use eirmod_calstr_buffered, only: deallocate_calstr_buffer
    deallocate(t_particles_loc)
    deallocate(t_calstr_loc)
    deallocate(t_postproc_loc)
    deallocate(n_particles_loc)
    deallocate(n_calstr_loc)
    if (my_pe == 0) then
      deallocate(t_particles_prev)
      deallocate(t_calstr_prev)
      deallocate(t_postproc_prev)
      if(allocated(throughput2)) then
        deallocate(throughput2)
        deallocate(t_overhead2)
      endif
    endif
    call deallocate_calstr_buffer
    return
  end subroutine deallocate_improved_strategy

  subroutine init_improved_strategy(ierror)
  ! Check if we have correct MPI version, and initializes calstr_buffer.
  ! This subroutine does not define a parallelization strategy. After calling
  ! this subroutine, please use any other strategy (preferably STRATEGY_APCAS)
  ! to initialize the work distribution.
    use eirmod_calstr_buffered
    use eirmod_mpi
    implicit none
    integer, intent(out) :: ierror !< 0 = success, any other values = error
    integer ierr, ierrr
    integer mpi_major, mpi_minor
    external eirene_masage, eirene_exit_own
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external mpi_allreduce
#endif

    ierror = 0

    ! At least one stratum is assumed to have particles.
    if (all(npts(1:nstrai) == 0)) then
      write(iunout,*) 'Error: no stratum with particles'
      ierror = 1
      return
    endif

    call mpi_get_version(mpi_major, mpi_minor, ierr)
    if (mpi_major < 3) then
      write(iunout,*) 'Error: balanced strategy requires MPI 3 '// &
                    & 'for non-blocking communication'
      ierror = 1
      return
    endif

    ! We cannot access my_pe and nprs from eirmod_cpes, because that would be a
    ! a circular reference. We store locally instead.
    call mpi_comm_rank(MPI_COMM_WORLD, my_pe, ierr)
    if (ierr /= mpi_success) then
      call eirene_masage('Error in init_improved_strategy at mpi_comm_rank.')
      call eirene_exit_own(1)
    end if
    call mpi_comm_size(MPI_COMM_WORLD, nprs, ierr)
    if (nprs == 1) then
      if (my_pe == 0) then
        write(iunout,*) 'Warning: improved strategy not needed in serial mode'
      endif
      ierror = 1
      return
    endif

    if (my_pe == 0) then
      write(iunout,*) 'Initializing improved strategy'
    endif
    call allocate_calstr_buffer(ierr)
    call mpi_allreduce(ierr, ierror, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierrr)
    if (ierror /= 0) then
      if(my_pe == 0) then
        write(iunout,*) 'Error: not enough memory to allocate calstr buffer'
      endif
    endif
    return
  end subroutine init_improved_strategy

  subroutine opt_improved_strategy(nparts_loc, npestr, stratum_leader, &
              procforstra)

    use eirmod_mpi
    implicit none
    integer, intent(out), dimension(:) :: nparts_loc
    integer, intent(out), dimension(:) :: npestr
    integer, intent(out), dimension(:) :: stratum_leader
    logical, intent(out), dimension(:,0:) :: procforstra

    ! mpi communication variables for performance counters
    integer :: n_calstr_sum(nstrai), n_particles_sum(nstrai), n_postproc_sum
    real(kind=dp), dimension(nstrai) :: t_calstr_sum, t_particles_sum, t_postproc_sum, t_calstr_ave
    real(kind=dp), dimension(0:nprs-1) :: t_calstr_pe, t_particles_pe, t_postproc_pe

    real(kind=dp), dimension(nstrai) :: throughput_strata
    integer, dimension(nstrai*nprs) :: nparts_loc_all
    real(kind=dp), dimension(0:nprs-1) :: t_pe
    real(kind=dp) :: t_average
    real(kind=dp) :: time
    real(kind=dp) :: dt
    integer, dimension(nstrai) :: npts_remaining
    integer :: i, j, k, idx, ierr, dn
    real(kind=dp) :: timetable(nstrai)
    integer :: mapping(nstrai), rank
    real(kind=dp) :: threshold
    ! reduction to avoid wait time of leader increasing with iterations
    real(kind=dp), parameter :: reduction = 0.75_dp
    integer, save :: iteration = 0

    external :: eirene_masage, eirene_exit_own
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external :: mpi_bcast, mpi_gather, mpi_reduce, mpi_scatter
#endif

    if (nprs == 1) return ! nothing to optimize for serial mode

    if (my_pe == 0) then
      write(iunout,*) 'Optimizing workload distribution'
    endif

    ! obtain the performace counters
    call mpi_reduce(n_calstr_loc, n_calstr_sum, nstrai, &
                    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(t_calstr_loc, t_calstr_sum, nstrai, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(n_particles_loc, n_particles_sum, nstrai, &
                    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(t_particles_loc, t_particles_sum, nstrai, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(n_postproc_loc, n_postproc_sum, 1, &
                    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(t_postproc_loc, t_postproc_sum, nstrai, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call mpi_gather(sum(t_calstr_loc), 1, MPI_DOUBLE_PRECISION, &
                    t_calstr_pe, 1, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    call mpi_gather(sum(t_particles_loc), 1, MPI_DOUBLE_PRECISION, &
                    t_particles_pe, 1, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    call mpi_gather(sum(t_postproc_loc), 1, MPI_DOUBLE_PRECISION, &
                    t_postproc_pe, 1, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    ! reset the counters to measure each iteration
    call reset_counters()

    if (my_pe == 0) then
      t_pe = 0.0
      npts_remaining = npts(1:nstrai)
      nparts_loc_all = 0
      stratum_leader = 0
      procforstra = .false.

      if (any(npts(1:nstrai) == 0)) then
        iteration = 0
      endif
      ! averaging for stability
      iteration = iteration + 1
      if (iteration == 1) then
        ! The first iteration uses results of APCAS.
        t_calstr_sum = 0.0_dp
        t_postproc_sum = sum(t_postproc_sum) / real(nstrai,dp)
      else if (iteration == 2) then
        ! Values of t_particles_sum is not stable yet.
        t_calstr_sum = 0.0_dp
        t_postproc_sum = (t_postproc_sum * 2.0_dp + t_postproc_prev) / 3.0_dp
      else if (iteration == 3) then
        t_calstr_sum = (t_calstr_sum * 2.0_dp + t_calstr_prev) / 3.0_dp
        t_postproc_sum = (t_postproc_sum * 2.0_dp + t_postproc_prev) / 3.0_dp
        t_particles_sum = (t_particles_sum * 2.0_dp + t_particles_prev) / 3.0_dp
      else
        t_calstr_sum = (t_calstr_sum + t_calstr_prev) / 2.0_dp
        t_postproc_sum = (t_postproc_sum + t_postproc_prev) / 2.0_dp
        t_particles_sum = (t_particles_sum + t_particles_prev) / 2.0_dp
      endif
      t_calstr_prev = t_calstr_sum
      t_postproc_prev = t_postproc_sum
      t_particles_prev = t_particles_sum

      do k = 1, nstrai
        ! averaged throughtput [particles/second] for each stratum
        if(t_particles_sum(k) == 0) then
          throughput_strata(k) = real(sum(n_particles_sum),dp) / sum(t_particles_sum)
        else
          throughput_strata(k) = real(n_particles_sum(k),dp) / t_particles_sum(k)
        endif
        ! averaged communication time within each stratum as a threshold time
        ! to decide separation of a stratum to the next process
        if(n_calstr_sum(k) == 0) then
          t_calstr_ave(k) = 0.0_dp
        else
          t_calstr_ave(k) = t_calstr_sum(k) / real(n_calstr_sum(k),dp) * reduction
        endif
      enddo
      if (sum(n_calstr_sum) == 0) then
        threshold = 0.0_dp
      else
        threshold = sum(t_calstr_sum) / real(sum(n_calstr_sum),dp)
      endif
      if (print_debug_info) then
        write(iunout,*) 'threshold time', threshold
      endif

      if (print_debug_info) then
        do k = 1, nstrai
          write(iunout,'(i4,i8,f8.3,i8,f8.3,f8.3,f8.3)') &
             &  k, n_particles_sum(k), t_particles_sum(k), n_calstr_sum(k), &
             &     t_calstr_sum(k), t_calstr_ave(k), t_postproc_sum(k)
        enddo
        do i = 0, nprs-1
          write(iunout,'(i4,f8.3,f8.3,f8.3)') &
             &  i, t_particles_pe(i), t_calstr_pe(i), t_postproc_pe(i)
        enddo
      endif  ! my_pe == 0

      ! averaged calculation time for each stratum
      timetable = t_particles_sum + t_postproc_sum + t_calstr_sum * reduction
      ! sort them in descending order
      call sort(nstrai, timetable, mapping)

      ! average calculation time per process
      t_average = (sum(t_particles_sum) + sum(t_postproc_sum)) / nprs
      do k = 1, nstrai
        time = t_particles_sum(k) + t_postproc_sum(k)
        if (time > t_average) then
          time = t_calstr_ave(k) * (floor(time / t_average) + 1.0_dp)
        else
          time = 0.0_dp
        endif
        t_average = t_average + time / nprs
      enddo
      if (threshold > (t_average / 4.0_dp) ) then
        write(iunout,*) &
          &  'Warning: too large threshold, reduced to t_average/4', &
          &   threshold, t_average
        threshold = t_average / 4.0_dp
      endif

      rank = 0
      ! allocating strata with no particles at first
      do j=1,nstrai
        if (npts_remaining(j) == 0) then
          procforstra(j,rank) = .true.
          stratum_leader(j) = rank
        endif
      enddo
      ! allocating strata with non-zero particles
      do j=1,nstrai
        ! average calculation time per process
        if (iteration /= 1) then
          t_average = 0.0_dp
          do k = 1, nstrai
            if (npts_remaining(k) /= 0) then
              t_average = t_average + t_particles_sum(k) + &
                        & t_postproc_sum(k) + t_calstr_sum(k) * reduction
            endif
          enddo
          t_average = (t_pe(rank) + t_average) / real(nprs - rank,dp)
        endif

        ! 1) choose the next stratum
        if (rank == 0) then
          ! processing the strata from the smallest one to keep many leaders on rank 0
          k = mapping(nstrai - j + 1)
          if (npts_remaining(k) == 0) then
            cycle
          endif
        else
          ! processing the longest stratum within the remaining time of each rank
          k = 0
          do idx=1,nstrai
            if (npts_remaining(mapping(idx)) == 0) then
              cycle
            endif
            k = mapping(idx)
            if (timetable(idx) < (t_average + threshold - t_pe(rank)) ) then
              exit
            endif
          enddo
          if (k == 0) then
            call eirene_masage('Error in opt_improved_strategy: k==0.')
            call eirene_exit_own(1)
          end if
        endif
        if (print_debug_info) then
          write(iunout,*) 'stratum',k
        endif

        ! 2) determine the leader and calculate the time per process for the strata
        time = t_pe(rank) + t_particles_sum(k) + t_postproc_sum(k) + t_calstr_sum(k)
        if (print_debug_info) then
          write(iunout,'(a,i4,i4,f8.3,f8.3)') 'k rank t_average time:', &
            k,rank,t_average,time
          write(iunout,'(a,f8.3,f8.3,f8.3)') 't_pe t_postproc_sum threshold:', &
            t_pe(rank),t_postproc_sum(k),threshold
        endif
        if (time < (t_average + threshold) ) then
          stratum_leader(k) = rank
          time = t_average + threshold
          if (print_debug_info) then
            write(iunout,'(a,i4)') 'stratum_leader(1)',stratum_leader(k)
          endif
        else
          if (rank == 0) then
            if (t_postproc_sum(k) < (t_average - t_pe(rank)) ) then
              stratum_leader(k) = rank
            else
              stratum_leader(k) = rank + 1
            endif
            if (print_debug_info) then
              write(iunout,'(a,i4)') 'stratum_leader(2)',stratum_leader(k)
            endif
          else if (t_pe(rank) == 0.0 .or. rank == nprs - 1) then
            stratum_leader(k) = rank
            if (print_debug_info) then
              write(iunout,'(a,i4,f8.3,f8.3,f8.3)') 'stratum_leader(3)', &
                 &  stratum_leader(k),t_pe(rank),time,t_average+threshold
            endif
          else
            stratum_leader(k) = rank + 1
            if (print_debug_info) then
              write(iunout,'(a,i4,i4,f8.3,f8.3)') 'stratum_leader(4)', &
                 &  stratum_leader(k)
            endif
          endif
          if (mod(time,t_average) < threshold) then
            time = time / floor(time / t_average)
          else if (t_average - mod(time,t_average) < threshold) then
            time = time / (floor(time / t_average) + 1.0_dp)
          else
            time = t_average
          endif
        endif
        t_pe(stratum_leader(k)) = t_pe(stratum_leader(k)) + t_postproc_sum(k)

        ! 3) assign the particles to MPI processes
        do while (npts_remaining(k) > 0)
          dt = time - t_calstr_ave(k) - t_pe(rank)
          dn = nint(dt * throughput_strata(k))
          if (dn <= 0) then
            if (stratum_leader(k) > rank) then
              rank = rank + 1
              cycle
            endif
            dn = 1
            dt = dn / throughput_strata(k)
          ! dn+1 is used instead of dn to adjust a numerical error
          else if (npts_remaining(k) <= dn + 1 .or. (j == nstrai .and. rank == nprs - 1)) then
            dn = npts_remaining(k)
            dt = dn / throughput_strata(k)
            if (stratum_leader(k) > rank) then
              dn = dn - 1
              dt = dn / throughput_strata(k)
            endif
          endif
          if (print_debug_info) then
            write(iunout,'(i4,i4,i8,f8.3,f8.3,f8.3,f8.3)') &
               &  k,rank,dn,dt,time,t_pe(rank),t_average
          endif
          idx = rank * nstrai + k
          t_pe(rank) = t_pe(rank) + dt
          if (dn /= npts(k)) then
            t_pe(rank) = t_pe(rank) + t_calstr_ave(k)
          endif
          nparts_loc_all(idx) = nparts_loc_all(idx) + dn
          npts_remaining(k) = npts_remaining(k) - dn
          procforstra(k,rank) = .true.
          if (print_debug_info) then
            write(iunout,'(a,i8,f8.3)') 'dn, dt =',dn,dt
          endif
          if ((npts_remaining(k) == 0 .and. t_pe(rank) > (t_average - threshold)) &
            .or. npts_remaining(k) /= 0) then
            if (rank == nprs - 1) then
              write(iunout,*) &
                 & 'Warning: particle assignment overflow '// &
                 & 'on the last rank at stratum', k, npts_remaining(k)
            else
              if (iteration == 1) then
                t_average = t_average + (t_average - t_pe(rank)) / (nprs - rank - 1)
              endif
              rank = rank + 1
            endif
          endif
        end do
      end do

      if (print_debug_info) then
        write(iunout,'(a,f8.3)') 'total time:',sum(t_pe)
        write(iunout,'(a,f8.3)') 'average time:',sum(t_pe)/nprs
        write(iunout,*) 'allocated time to each rank:'
        do rank = 0, nprs-1
          write(iunout,'(i4,f8.3,f8.3)') rank, t_pe(rank), t_pe(rank)-sum(t_pe)/nprs
        enddo
      endif

      ! 4) calculate the overhead for output message and the number or processes for each stratum
      if (.not. allocated(throughput2)) then
        allocate(throughput2(nstrai))
        allocate(t_overhead2(0:nprs-1))
      endif
      throughput2 = throughput_strata
      t_overhead2 = 0.0_dp
      do k = 1, nstrai
        rank = stratum_leader(k)
        t_overhead2(rank) = t_overhead2(rank) + t_postproc_sum(k)
        npestr(k) = 0
        do rank = 0, nprs-1
          idx = rank * nstrai + k
          if (nparts_loc_all(idx) > 0) then
            npestr(k) = npestr(k) + 1
            t_pe(rank) = t_pe(rank) - t_calstr_ave(k)
            t_overhead2(rank) = t_overhead2(rank) + t_calstr_ave(k)
          endif
        end do
      end do

    endif  ! my_pe == 0

    call mpi_bcast(stratum_leader, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npestr, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_scatter(nparts_loc_all, nstrai, MPI_INTEGER, &
                 &   nparts_loc, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! procforstra(nstra,0:nprs-1): nstra may different from nstrai.
    call mpi_bcast(procforstra, size(procforstra), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    return
  end subroutine opt_improved_strategy

  subroutine reset_counters()
    n_calstr_loc = 0
    t_calstr_loc = 0.0_dp
    n_particles_loc = 0
    t_particles_loc = 0.0_dp
    n_postproc_loc = 0
    t_postproc_loc = 0.0_dp
    return
  end subroutine reset_counters

  subroutine time_calstr_improved(istra, time)
    integer, intent(in) :: istra
    real(kind=dp), intent(in) :: time

    if(time /= 0.0_dp) then
      t_calstr_loc(istra) = t_calstr_loc(istra) + time
      n_calstr_loc(istra) = n_calstr_loc(istra) + 1
    endif
    return
  end subroutine time_calstr_improved

  subroutine time_particles_improved(istra, time, n)
    integer, intent(in) :: istra !< stratum idx
    real(kind=dp), intent(in) :: time !< time (s)
    integer, intent(in) :: n !< number of particles

    t_particles_loc(istra) = t_particles_loc(istra) + time
    n_particles_loc(istra) = n_particles_loc(istra) + n
    return
  end subroutine time_particles_improved

  subroutine time_postproc_improved(istra, time)
    integer, intent(in) :: istra !< stratum idx
    real(kind=dp), intent(in) :: time !< time (s)

    t_postproc_loc(istra) = t_postproc_loc(istra) + time
    n_postproc_loc = n_postproc_loc + 1
    return
  end subroutine time_postproc_improved

  ! sort the elements in the array in descending order and make an index mapping
  subroutine sort(size, value, index)
    implicit none
    integer, intent(in) :: size
    real(kind=dp), intent(inout) :: value(1:size)
    integer, intent(out) :: index(1:size)
    integer :: first, last, i, j, itmp
    real(kind=dp) :: rtmp

    do i = 1, size
      index(i) = i
    enddo
    first = 1
    last = size
    do while(first < last)
      j = first
      do i = first, last - 1
        if(value(i) < value(i+1)) then
          rtmp = value(i)
          value(i) = value(i+1)
          value(i+1) = rtmp
          itmp = index(i)
          index(i) = index(i+1)
          index(i+1) = itmp
          j = i
        endif
      enddo
      last = j
    enddo
    return
  end subroutine sort

end module eirmod_improved_strategy
