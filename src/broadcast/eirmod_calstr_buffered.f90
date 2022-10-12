module eirmod_calstr_buffered
! A buffered version of calstr to allow non-blocking reduction
!
! All the data that is communicated in calstr are packed into a few large
! buffers (one for each data type), and reduced using MPI_IREDUCE.
!
! If you add any variables to calstr, then the same has to be added to this
! module too. Please add them to the following three subroutines:
!  - allocate_calstr_buffer,
!  - calstr_pack_all,
!  - calstr_unpack_all
!
! Note that the reduction operation is not progressing magically in the
! background after we call MPI_IREDUCE. We have to check the status of the
! reductions regularily using MPI_TEST, and this ensures that the reduction is
! progressing. This is implemented in calstr_progress_message.
!
! Alternatively one could use a thread in the background to manage the MPI
! communication (Intel MPI does that automatically if we set
! I_MPI_ASYNC_PROGRESS=1, but that was 40% slower using intel MPI 5.1).
!
! This module has one knob for tuning the communication: n_progress_check.
! This tells how many times we call MPI_TEST during the particle loop.

  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_cestim
  use eirmod_cspez
  use eirmod_comprt
  use eirmod_csdvi
  use eirmod_coutau
  use eirmod_mpi
  use eirmod_cstep
  use eirmod_cai

  implicit none

  private
  public :: allocate_calstr_buffer, deallocate_calstr_buffer
  public :: calstr_buffered_finish
  public :: calstr_progress_message
  public :: eirene_calstr_buffered

  !> nubmer of buffers (1 double and 1 logical)
  integer, parameter :: N_BUFFERS = 2
  !> Identifies the communication operations, one request object for each buffer
  integer, save, dimension(N_BUFFERS) :: calstr_request =  MPI_REQUEST_NULL

  real(kind=dp), allocatable, dimension(:) :: calstr_buffer_d
  integer :: pos_d !< position in the double precision buffer

  logical, allocatable, dimension(:) :: calstr_buffer_l
  integer :: pos_l !< position in the logical buffer

  ! To ensure that the non-blocking reduce operations are progressing, we
  ! periodically call MPI_TEST for the calstr_request array.
  ! We would like to perform n_progress_check test.
  integer, parameter :: n_progress_check = 100

  ! For each stratum we process 1..nparst_loc(istra) particles
  ! so we will test after every nparts_loc(istra)/n_progress_check particles.
  ! Next_check(istra) stores the particle number when the next check is due
  integer, allocatable, dimension(:) :: next_check

  contains

  subroutine calstr_init_progress
    use eirmod_mpi
    integer :: ierr
    if(.not.allocated (next_check)) then
      write(iunout,*) 'Error next_check should be allocated'
      call mpi_abort(mpi_comm_world, -1, ierr)
    endif
    next_check = 0
  end subroutine calstr_init_progress

  subroutine calstr_progress_message(n, istra, nparts_loc)
  ! Performs an MPI_TEST on the non-blocking messages, to ensure that they are
  ! progressing.
  !
  ! The calstr reduction for stratum(k-1) should finish before we start calstr
  ! for stratum(k). Therefore, during the particle loop of stratum(k) we
  ! periodically call MPI_TEST. We aim to call the test subroutine
  ! n_progress_check times. Using the particle loop index n, and the total
  ! number of particles nparts_loc(k) we decide whether to call MPI_TEST or
  ! skip it.
    use eirmod_mpi
    integer, intent(in) :: n !< particle index
    integer, intent(in) :: istra !< stratum number
    !> total number of particles
    integer, intent(in), dimension(:), allocatable :: nparts_loc
    logical :: flag
    integer :: ierr
    if (.not. allocated(next_check)) then
      ! We are not in buffered mode in this case
      return
    endif
    if (n > next_check(istra)) then
      call MPI_Testall(N_BUFFERS, calstr_request, flag, MPI_STATUSES_IGNORE, ierr)
      next_check(istra) = next_check(istra) + &
                          max(1, nparts_loc(istra)/n_progress_check)
    endif
  end subroutine calstr_progress_message

  subroutine calstr_buffered_finish()
    use eirmod_mpi
    integer :: ierr
    call MPI_WAITALL(N_BUFFERS, calstr_request, MPI_STATUSES_IGNORE, ierr)
  end subroutine calstr_buffered_finish

  subroutine eirene_calstr_buffered(comm, my_pe_gr, leader)
  ! Non-blocking reduction using extra buffers.
  ! The data is first copied to buffers using calstr_pack_all, and the buffers
  ! are reduced using MPI_IREDUCE. The stratum leader waits for the operation
  ! to finish, the others are free to continue with the next stratum.
  !
  ! If you need to add any other reduction operatons, then please do so in
  ! allocate_calstr_buffer, calstr_pack_all, calstr_unpack_all.
    use eirmod_mpi
    integer, intent(in) :: comm !< calstr_communicator for the stratum
    integer, intent(in) :: my_pe_gr !< my rank within the calstr communicator
    logical, intent(in) :: leader !< whether the PE is the leader of the stratum
    integer :: ierr, isize
    call calstr_init_progress
    ! Wait if there are any pending reduce operations for the previous stratum
    call MPI_WAITALL(N_BUFFERS, calstr_request, MPI_STATUSES_IGNORE, ierr)
    ! The previous calls have finished, the buffers are ready to use
    call calstr_pack_all
    if (leader) then
      if (my_pe_gr.ne.0) then
        ! During the creation of the stratum communicators we enforce that the
        ! leader will have rank 0. But it cannot hurt to cross-check it.
        write (iunout,*) 'Error, the rank of the stratum leader ', &
                     'must be zero, within the stratum.'
        call mpi_comm_size(comm, isize, ierr)
        write (0,*) 'eirene_calstr_buffered: my_pe_gr, size = ', my_pe_gr, isize
        call mpi_abort(MPI_COMM_WORLD, -1, ierr)
      endif
      call mpi_ireduce(MPI_IN_PLACE, calstr_buffer_d, &
             size(calstr_buffer_d), MPI_DOUBLE_PRECISION, MPI_SUM, &
             0, comm, calstr_request(1), ierr)
      if (ierr /= mpi_success) then
        call eirene_masage('ERROR IN eirene_calstr_buffered for calstr_buffer_d')
        call eirene_exit_own(1)
      end if
      call mpi_ireduce(MPI_IN_PLACE, calstr_buffer_l, &
             size(calstr_buffer_l), MPI_LOGICAL, MPI_LOR, &
             0, comm, calstr_request(2), ierr)
      if (ierr /= mpi_success) then
        call eirene_masage('ERROR IN eirene_calstr_buffered for calstr_buffer_l')
        call eirene_exit_own(1)
      end if
      ! The stratum leader has to postprocess the results, so it waits for the
      ! communication to finish
      call MPI_WAITALL(N_BUFFERS, calstr_request, MPI_STATUSES_IGNORE, ierr)
      call calstr_unpack_all
    else
      call mpi_ireduce(calstr_buffer_d, calstr_buffer_d, &
             size(calstr_buffer_d), MPI_DOUBLE_PRECISION, MPI_SUM, &
             0, comm, calstr_request(1), ierr)
      call mpi_ireduce(calstr_buffer_l, calstr_buffer_l, &
             size(calstr_buffer_l), MPI_LOGICAL, MPI_LOR, &
             0, comm, calstr_request(2), ierr)
    endif
  end subroutine eirene_calstr_buffered

  subroutine calstr_pack_d0(buff)
    real(kind=dp), intent(in) :: buff
    if (.not. allocated(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer_d not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_d + 1 > size(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer_d size is too small', &
                  pos_d + 1,  size(calstr_buffer_d)
      call eirene_exit_own(1)
    endif
    calstr_buffer_d(pos_d+1) = buff
    pos_d = pos_d + 1
  end subroutine calstr_pack_d0

  subroutine calstr_pack_d(n, buff)
    integer :: n
    real(kind=dp), dimension(*), intent(in) :: buff

    if (.not. allocated(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer_d not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_d + n > size(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer_d size is too small', &
                  pos_d + n,  size(calstr_buffer_d)
      call eirene_exit_own(1)
    endif
    calstr_buffer_d(pos_d+1:pos_d+n) = buff(1:n)
    pos_d = pos_d + n
  end subroutine calstr_pack_d

  subroutine calstr_unpack_d0(buff)
    real(kind=dp), intent(out) :: buff
    if (.not. allocated(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_d + 1 > size(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer size is too small', &
                  pos_d + 1,  size(calstr_buffer_d)
      call eirene_exit_own(1)
    endif
    buff = calstr_buffer_d(pos_d+1)
    pos_d = pos_d + 1
  end subroutine calstr_unpack_d0

  subroutine calstr_unpack_d(n, buff)
    integer, intent(in) :: n
    real(kind=dp), dimension(*) :: buff
    if (.not. allocated(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_d + n > size(calstr_buffer_d)) then
      write(iunout,*) 'Error, calstr_buffer size is too small', &
                  pos_d + n,  size(calstr_buffer_d)
      call eirene_exit_own(1)
    endif
    buff(1:n) = calstr_buffer_d(pos_d+1:pos_d+n)
    pos_d = pos_d + n
  end subroutine calstr_unpack_d

  subroutine calstr_pack_l0(buff)
    logical, intent(in) :: buff
    if (.not. allocated(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer_l not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_l + 1 > size(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer_l size is too small', &
                  pos_l + 1,  size(calstr_buffer_l)
      call eirene_exit_own(1)
    endif
    calstr_buffer_l(pos_l+1) = buff
    pos_l = pos_l + 1
  end subroutine calstr_pack_l0

  subroutine calstr_pack_l(n, buff)
    integer, intent(in) :: n
    logical, dimension(*), intent(in) :: buff

    if (.not. allocated(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer_l not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_l + n > size(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer_l size is too small', &
                  pos_l + n,  size(calstr_buffer_l)
      call eirene_exit_own(1)
    endif
    calstr_buffer_l(pos_l+1:pos_l+n) = buff(1:n)
    pos_l = pos_l + n
  end subroutine calstr_pack_l

  subroutine calstr_unpack_l0(buff)
    logical, intent(out) :: buff
    if (.not. allocated(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_l + 1 > size(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer size is too small', &
                  pos_l + 1,  size(calstr_buffer_l)
      call eirene_exit_own(1)
    endif
    buff = calstr_buffer_l(pos_l+1)
    pos_l = pos_l + 1
  end subroutine calstr_unpack_l0

  subroutine calstr_unpack_l(n, buff)
    integer, intent(in) :: n
    logical, dimension(*), intent(out) :: buff

    if (.not. allocated(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer not allocated'
      call eirene_exit_own(1)
    endif
    if (pos_l + n > size(calstr_buffer_l)) then
      write(iunout,*) 'Error, calstr_buffer size is too small', &
                  pos_l + n,  size(calstr_buffer_l)
      call eirene_exit_own(1)
    endif
    buff(1:n) = calstr_buffer_l(pos_l+1:pos_l+n)
    pos_l = pos_l + n
  end subroutine

  subroutine allocate_calstr_buffer(ierr)
    integer, intent(out) :: ierr
    integer n, ns, ispc
    integer :: istra
    istra = 1
    n = 0
    n = n + size(WTOTM(0:nmoli,istra))
    n = n + size(WTOTA(0:natmi,istra))
    n = n + size(WTOTI(0:nioni,istra))
    n = n + size(WTOTP(0:nplsi,istra))
    n = n + 1 ! = size(WTOTE(istra)), but size only works for arrays
    n = n + 1 ! size(XMCP(istra))
    n = n + 1 ! size(XMCT(istra))
    n = n + 1 ! size(PTRASH(istra))
    n = n + 1 ! size(ETRASH(istra))
    n = n + 1 ! size(ETOTA(istra))
    n = n + 1 ! size(ETOTM(istra))
    n = n + 1 ! size(ETOTI(istra))
    n = n + 1 ! size(ETOTP(istra))
    n = n + size(EELFI(0:nioni,istra))
!  particle balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    n = n + size(PPATI(0:natmi,istra))
    n = n + size(PPMLI(0:nmoli,istra))
    n = n + size(PPIOI(0:nioni,istra))
    if (nphoti > 0) then
      n = n + size(PPPHTI(0:nphoti,istra))
    end if
    n = n + size(PPPLI(0:nplsi,istra))
!  energy balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    n = n + 1 ! size(EPATI(istra))
    n = n + 1 ! size(EPMLI(istra))
    n = n + 1 ! size(EPIOI(istra))
    n = n + 1 ! size(EPPHTI(istra))
    n = n + size(EPPLI(0:nplsi,istra))
!  all other volume-averaged tallies: estimv
    n = n + size(ESTIMV)
!  all surface-averaged tallies: estims
    n = n + size(ESTIMS)
!   energy-resolved ("spectra") tallies
    do ispc=1,nadspc
      ns = estiml(ispc)%nspc
      n = n + ns + 2
!  standard deviation of energy-resolved "spectra"
      if (nsigi_spc > 0) then
       n = n + ns + 2
       n = n + ns + 2
       n = n + 1 ! size(estiml(ispc)%sgms)
      end if
    end do
!  standard deviation of volume-averaged tallies
    if (nsd > 0) then
      n = n + size(SDVI1)
    end if
!  standard deviation of surface-averaged tallies
    if (nsdw > 0) then
      n = n + size(SDVI2)
    end if
!  covariances between two volume-averaged tallies
    if (ncv > 0) then
       n = n + size(SIGMAC)
       n = n + size(SGMCS)
    end if

    allocate(next_check(nstra))
    allocate(calstr_buffer_d(n), stat=ierr)
    if (ierr.ne.0) then
      write(iunout,*) &
        'Error, not enough memory to allocate calstr_buffer_d with ', &
         n, ' elements'
      return
    end if
! logical arrays
    n = 0
    n = n + size(LOGMOL(0:nmoli,ISTRA))
    n = n + size(LOGATM(0:natmi,ISTRA))
    n = n + size(LOGION(0:nioni,ISTRA))
    if (nphoti > 0) then
      n = n + size(LOGPHOT(0:nphoti,ISTRA))
    end if
    n = n + size(LOGPLS(0:nplsi,ISTRA))

    allocate(calstr_buffer_l(n), stat=ierr)
    if (ierr.ne.0) then
      write(iunout,*) &
        'Error, not enough memory to allocate calstr_buffer_l with ', &
         n, ' elements'
      return
    end if
    ! variables from call eirene_calstr_usr should be also included, if there are any
  end subroutine allocate_calstr_buffer

  subroutine deallocate_calstr_buffer
    if (allocated(calstr_buffer_d)) deallocate(calstr_buffer_d)
    if (allocated(calstr_buffer_l)) deallocate(calstr_buffer_l)
    if (allocated(next_check)) deallocate(next_check)
  end subroutine deallocate_calstr_buffer

  subroutine calstr_pack_all
    use eirmod_mpi
    integer :: ns, ispc
    if (any(calstr_request.ne.MPI_REQUEST_NULL)) then
      write(iunout,*) 'Error, buffers are not ready for use'
      call MPI_ABORT(MPI_COMM_WORLD, -1, ispc)
    endif
    pos_d=0
    call calstr_pack_d(size(WTOTM(0:nmoli,istra)), WTOTM(0:nmoli,istra))
    call calstr_pack_d(size(WTOTA(0:natmi,istra)), WTOTA(0:natmi,istra))
    call calstr_pack_d(size(WTOTI(0:nioni,istra)), WTOTI(0:nioni,istra))
    call calstr_pack_d(size(WTOTP(0:nplsi,istra)), WTOTP(0:nplsi,istra))
    call calstr_pack_d(1, WTOTE(istra))
    call calstr_pack_d(1, XMCP(istra))
    call calstr_pack_d(1, XMCT(istra))
    call calstr_pack_d(1, PTRASH(istra))
    call calstr_pack_d(1, ETRASH(istra))
    call calstr_pack_d(1, ETOTA(istra))
    call calstr_pack_d(1, ETOTM(istra))
    call calstr_pack_d(1, ETOTI(istra))
    call calstr_pack_d(1, ETOTP(istra))
    call calstr_pack_d(size(EELFI(0:nioni,istra)), EELFI(0:nioni,istra))
!  particle balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    call calstr_pack_d(size(PPATI(0:natmi,istra)), PPATI(0:natmi,istra))
    call calstr_pack_d(size(PPMLI(0:nmoli,istra)), PPMLI(0:nmoli,istra))
    call calstr_pack_d(size(PPIOI(0:nioni,istra)), PPIOI(0:nioni,istra))
    if (nphoti > 0) then
      call calstr_pack_d(size(PPPHTI(0:nphoti,istra)), PPPHTI(0:nphoti,istra))
    end if
    call calstr_pack_d(size(PPPLI(0:nplsi,istra)), PPPLI(0:nplsi,istra))
!  energy balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    call calstr_pack_d(1, EPATI(istra))
    call calstr_pack_d(1, EPMLI(istra))
    call calstr_pack_d(1, EPIOI(istra))
    call calstr_pack_d(1, EPPHTI(istra))
    call calstr_pack_d(size(EPPLI(0:nplsi,istra)), EPPLI(0:nplsi,istra))
!  all other volume-averaged tallies: estimv
    call calstr_pack_d(size(ESTIMV), ESTIMV)
!  all surface-averaged tallies: estims
    call calstr_pack_d(size(ESTIMS), ESTIMS)
!   energy-resolved ("spectra") tallies
    do ispc=1,nadspc
      ns = estiml(ispc)%nspc
      call calstr_pack_d(ns+2, estiml(ispc)%spc(0:ns+1))
!  standard deviation of energy-resolved "spectra"
      if (nsigi_spc > 0) then
        call calstr_pack_d(ns+2, estiml(ispc)%sdv(0:ns+1))
        call calstr_pack_d(ns+2, estiml(ispc)%sgm(0:ns+1))
        call calstr_pack_d0(estiml(ispc)%sgms)
      end if
    end do
!  standard deviation of volume-averaged tallies
    if (nsd > 0) then
      call calstr_pack_d(size(SDVI1), SDVI1)
    end if
!  standard deviation of surface-averaged tallies
    if (nsdw > 0) then
      call calstr_pack_d(size(SDVI2), SDVI2)
    end if
!  covariances between two volume-averaged tallies
    if (ncv > 0) then
      call calstr_pack_d(size(SIGMAC), SIGMAC)
      call calstr_pack_d(size(SGMCS), SGMCS)
    end if

    pos_l=0
    ! logical arrays
    call calstr_pack_l(size(LOGMOL(0:nmoli,ISTRA)), LOGMOL(0:nmoli,ISTRA))
    call calstr_pack_l(size(LOGATM(0:natmi,ISTRA)), LOGATM(0:natmi,ISTRA))
    call calstr_pack_l(size(LOGION(0:nioni,ISTRA)), LOGION(0:nioni,ISTRA))
    if (nphoti > 0) then
      call calstr_pack_l(size(LOGPHOT(0:nphoti,ISTRA)), LOGPHOT(0:nphoti,ISTRA))
    end if
    call calstr_pack_l(size(LOGPLS(0:nplsi,ISTRA)), LOGPLS(0:nplsi,ISTRA))
    ! variables from eirene_calstr_usr should be also packed
  end subroutine calstr_pack_all

  subroutine calstr_unpack_all
    use eirmod_mpi
    integer :: ns, ispc
    if (any(calstr_request.ne.MPI_REQUEST_NULL)) then
      write(iunout,*) 'Error, buffers are not ready for use'
      call MPI_ABORT(MPI_COMM_WORLD, -1, ispc)
    endif
    pos_d=0
    call calstr_unpack_d(size(WTOTM(0:nmoli,istra)), WTOTM(0:nmoli,istra))
    call calstr_unpack_d(size(WTOTA(0:natmi,istra)), WTOTA(0:natmi,istra))
    call calstr_unpack_d(size(WTOTI(0:nioni,istra)), WTOTI(0:nioni,istra))
    call calstr_unpack_d(size(WTOTP(0:nplsi,istra)), WTOTP(0:nplsi,istra))
    call calstr_unpack_d(1, WTOTE(istra))
    call calstr_unpack_d(1, XMCP(istra))
    call calstr_unpack_d(1, XMCT(istra))
    call calstr_unpack_d(1, PTRASH(istra))
    call calstr_unpack_d(1, ETRASH(istra))
    call calstr_unpack_d(1, ETOTA(istra))
    call calstr_unpack_d(1, ETOTM(istra))
    call calstr_unpack_d(1, ETOTI(istra))
    call calstr_unpack_d(1, ETOTP(istra))
    call calstr_unpack_d(size(EELFI(0:nioni,istra)), EELFI(0:nioni,istra))
!  particle balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    call calstr_unpack_d(size(PPATI(0:natmi,istra)), PPATI(0:natmi,istra))
    call calstr_unpack_d(size(PPMLI(0:nmoli,istra)), PPMLI(0:nmoli,istra))
    call calstr_unpack_d(size(PPIOI(0:nioni,istra)), PPIOI(0:nioni,istra))
    if (nphoti > 0) then
      call calstr_unpack_d(size(PPPHTI(0:nphoti,istra)), PPPHTI(0:nphoti,istra))
    end if
    call calstr_unpack_d(size(PPPLI(0:nplsi,istra)), PPPLI(0:nplsi,istra))
!  energy balance tallies: from bulk (ipls) to species a,m,i,ph,pl
    call calstr_unpack_d(1, EPATI(istra))
    call calstr_unpack_d(1, EPMLI(istra))
    call calstr_unpack_d(1, EPIOI(istra))
    call calstr_unpack_d(1, EPPHTI(istra))
    call calstr_unpack_d(size(EPPLI(0:nplsi,istra)), EPPLI(0:nplsi,istra))
!  all other volume-averaged tallies: estimv
    call calstr_unpack_d(size(ESTIMV), ESTIMV)
!  all surface-averaged tallies: estims
    call calstr_unpack_d(size(ESTIMS), ESTIMS)
!   energy-resolved ("spectra") tallies
    do ispc=1,nadspc
      ns = estiml(ispc)%nspc
      call calstr_unpack_d(ns+2, estiml(ispc)%spc(0:ns+1))
!  standard deviation of energy-resolved "spectra"
      if (nsigi_spc > 0) then
        call calstr_unpack_d(ns+2, estiml(ispc)%sdv(0:ns+1))
        call calstr_unpack_d(ns+2, estiml(ispc)%sgm(0:ns+1))
        call calstr_unpack_d0(estiml(ispc)%sgms)
      end if
    end do
!  standard deviation of volume-averaged tallies
    if (nsd > 0) then
      call calstr_unpack_d(size(SDVI1), SDVI1)
    end if
!  standard deviation of surface-averaged tallies
    if (nsdw > 0) then
      call calstr_unpack_d(size(SDVI2), SDVI2)
    end if
!  covariances between two volume-averaged tallies
    if (ncv > 0) then
       call calstr_unpack_d(size(SIGMAC), SIGMAC)
       call calstr_unpack_d(size(SGMCS), SGMCS)
    end if
! logical arrays
    pos_l = 0
    call calstr_unpack_l(size(LOGMOL(0:nmoli,ISTRA)), LOGMOL(0:nmoli,ISTRA))
    call calstr_unpack_l(size(LOGATM(0:natmi,ISTRA)), LOGATM(0:natmi,ISTRA))
    call calstr_unpack_l(size(LOGION(0:nioni,ISTRA)), LOGION(0:nioni,ISTRA))
    if (nphoti > 0) then
      call calstr_unpack_l(size(LOGPHOT(0:nphoti,ISTRA)), LOGPHOT(0:nphoti,ISTRA))
    end if
    call calstr_unpack_l(size(LOGPLS(0:nplsi,ISTRA)), LOGPLS(0:nplsi,ISTRA))
    ! variables from eirene_calstr_usr should be also unpacekd
  end subroutine calstr_unpack_all


end module eirmod_calstr_buffered
