!> Wrapper module for MPI, and dummy module for serial compilation
module eirmod_mpi
#ifdef USE_MPI
  implicit none
  include 'mpif.h'

#if MPI_VERSION < 3
! MPI libraries with MPI version 3 are available on all platforms.
! If anyone still wants to use an older library without MPI 3 subroutines,
! then here we need to define a dummy MPI_IREDUCE
  interface mpi_ireduce
    ! During execution, we check the MPI version to avoid calling these if an
    ! older MPI library is used. If we still call MPI_IREDUCE, then it is an
    ! error. These dummy implementations write the error message and abort
    ! execution.
#ifndef GFORTRAN
    module procedure mpi_ireduce_i0_r1
    module procedure mpi_ireduce_r1_r1
    module procedure mpi_ireduce_i0_l1
    module procedure mpi_ireduce_l1_l1
#endif
  end interface

  contains
#endif

#else

  ! Dummy MPI module for serial compilation
  implicit none
  integer, parameter :: mpi_success = 0
  integer, parameter :: mpi_failure = 1
  integer, parameter :: mpi_comm_world = 0
  integer, parameter :: mpi_comm_self = 0
  integer, parameter :: mpi_comm_null = -1
  integer, parameter :: mpi_undefined = -32766

  integer, parameter :: mpi_group_null = mpi_comm_null
  integer, parameter :: mpi_group_empty = mpi_group_null

!  Dummy of an include file needed by MPI
  integer, parameter :: MPI_INTEGER = 1,   &
                 MPI_REAL = 2,             &
                 MPI_DOUBLE_PRECISION = 3, &
                 MPI_REAL8 = 3,            &
                 MPI_LOGICAL = 4,          &
                 MPI_CHARACTER = 5,        &
                 MPI_STATUS_IGNORE = -1,   &
                 MPI_STATUSES_IGNORE = -1, &
                 MPI_STATUS_SIZE = 3,      &
                 MPI_SUM = 1,              &
                 MPI_LOR = 2,              &
                 MPI_IN_PLACE = 0

  integer, parameter :: MPI_REQUEST_NULL = 0
  integer, parameter :: MPI_VERSION = 2
  integer, parameter :: MPI_SUBVERSION = 2


  interface mpi_allgather
    module procedure mpi_allgather_i0, mpi_allgather_i1
    module procedure mpi_allgather_l0
  end interface

  interface mpi_allreduce
    module procedure mpi_allreduce_i0
    module procedure mpi_allreduce_r0
  end interface

  interface mpi_gather
     module procedure mpi_gather_i0_i1
     module procedure mpi_gather_i1_i1
     module procedure mpi_gather_r0_r0
     module procedure mpi_gather_r0_r1
     module procedure mpi_gather_r1_r1
  end interface

  interface mpi_gatherv
     module procedure mpi_gatherv_i2_i2
     module procedure mpi_gatherv_r2_r2
  end interface

  interface mpi_scatter
     module procedure mpi_scatter_i0_i0
     module procedure mpi_scatter_i1_i0
     module procedure mpi_scatter_i1_i1
  end interface

  interface mpi_ireduce
    ! During execution, we check the MPI version to avoid calling these if an
    ! older MPI library is used. If we still call MPI_IREDUCE, then it is an
    ! error. These dummy implementations write the error message and abort
    ! execution.
    module procedure mpi_ireduce_i0_r1
    module procedure mpi_ireduce_r1_r1
    module procedure mpi_ireduce_i0_l1
    module procedure mpi_ireduce_l1_l1
  end interface

  interface mpi_send
    module procedure mpi_send_dum_a0, mpi_send_dum_a1, mpi_send_dum_a2
    module procedure mpi_send_dum_a3
    module procedure mpi_send_dum_i0, mpi_send_dum_i1, mpi_send_dum_i2
  end interface

  interface mpi_recv
    module procedure mpi_recv_dum_a0, mpi_recv_dum_a1, mpi_recv_dum_a2
    module procedure mpi_recv_dum_a3
    module procedure mpi_recv_dum_i0, mpi_recv_dum_i1, mpi_recv_dum_i2
  end interface

  interface mpi_reduce
    module procedure mpi_reduce_i0_i0, mpi_reduce_i1_i1
    module procedure mpi_reduce_L1_L1, mpi_reduce_L2_L1
    module procedure mpi_reduce_inplace_r0, mpi_reduce_inplace_r1, mpi_reduce_inplace_r2, &
      mpi_reduce_inplace_r3, mpi_reduce_inplace_r4
    module procedure mpi_reduce_r0_r0, mpi_reduce_r1_r1, mpi_reduce_r2_r2, &
      mpi_reduce_r3_r3, mpi_reduce_r4_r4
  end interface

  interface mpi_bcast
      module procedure MPI_BCAST_DUM_L, MPI_BCAST_DUM_L1, MPI_BCAST_DUM_L2
      module procedure MPI_BCAST_DUM_R_DP
      module procedure MPI_BCAST_DUM_A1, MPI_BCAST_DUM_A2
      module procedure MPI_BCAST_DUM_A3, MPI_BCAST_DUM_A4
      module procedure MPI_BCAST_DUM_A5, MPI_BCAST_DUM_A6
      module procedure MPI_BCAST_DUM_I0, MPI_BCAST_DUM_I1
      module procedure MPI_BCAST_DUM_I2, MPI_BCAST_DUM_I3
      module procedure MPI_BCAST_DUM_I6
      module procedure MPI_BCAST_DUM_C0, MPI_BCAST_DUM_C1
      module procedure MPI_BCAST_DUM_C2
  end interface

  interface mpi_testall
    module procedure mpi_testall_i0, mpi_testall_i1
  end interface

  interface mpi_waitall
    module procedure mpi_waitall_i0, mpi_waitall_i1
  end interface

  contains

  subroutine mpi_abort(comm, errorcode, ierr)
    use eirmod_comprt
    implicit none
    integer, intent(in) :: comm, errorcode
    integer, intent(out) :: ierr
    ierr = MPI_SUCCESS
    write (iunout,*) 'MPI ABORT: Shut down with error code: ', errorcode
    stop
  end subroutine

  subroutine mpi_bcast_dum_l (buffer,cnt,datatype,root,comm,ier)
    implicit none
    integer, intent(out) :: ier
    logical :: buffer
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_l1 (buffer,cnt,datatype,root,comm,ier)
    implicit none
    integer, intent(out) :: ier
    logical :: buffer(:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_l2 (buffer,cnt,datatype,root,comm,ier)
    implicit none
    integer, intent(out) :: ier
    logical :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_r_dp (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp)  :: buffer
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a1 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a2 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a3 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a4 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a5 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_a6 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:,:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_i0 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer  :: buffer
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_i1 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

    subroutine mpi_bcast_dum_i2 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer  :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_i3 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_i6 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:,:,:,:,:,:)
    integer, intent(in) :: cnt,datatype,root,comm
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_c0 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer, intent(in) :: cnt,datatype,root,comm
    character(*) :: buffer
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_c1 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer, intent(in) :: cnt,datatype,root,comm
    character(*) :: buffer(*)
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_bcast_dum_c2 (buffer,cnt,datatype,root,comm,ier)
    use eirmod_precision
    implicit none
    integer, intent(out) :: ier
    integer, intent(in) :: cnt,datatype,root,comm
    character(*) :: buffer(1,*)
    ier = MPI_SUCCESS
  end subroutine

  subroutine mpi_allgather_i0 (data1, nsend, sendtype, data2, nrecv, recvtype, &
    comm, ierror )
  integer, intent(in) :: nsend
  integer, intent(in) :: data1
  integer, intent(out) :: data2(:)
  integer, intent(in) :: comm, nrecv, recvtype, sendtype
  integer, intent(out) :: ierror
  ierror = MPI_SUCCESS
  data2(1) = data1
  end subroutine

  subroutine mpi_allgather_i1 (data1, nsend, sendtype, data2, nrecv, recvtype, &
    comm, ierror )
  integer, intent(in) :: nsend
  integer, intent(in) :: data1(:)
  integer, intent(out) :: data2(:)
  integer, intent(in) :: comm, nrecv, recvtype, sendtype
  integer, intent(out) :: ierror
  if (nsend==nrecv) then
    ierror = MPI_SUCCESS
    data2(1:nsend) = data1(1:nsend)
  else
    ierror = MPI_FAILURE
  endif
  end subroutine

  subroutine mpi_allgather_l0 (data1, nsend, sendtype, data2, nrecv, recvtype, &
    comm, ierror )
  integer, intent(in) :: nsend
  logical, intent(in) :: data1
  logical, intent(out) :: data2(:)
  integer, intent(in) :: comm, nrecv, recvtype, sendtype
  integer, intent(out) :: ierror
  ierror = MPI_SUCCESS
  data2(1) = data1
  end subroutine

  subroutine mpi_allreduce_i0 (data1, data2, n, dtype, op, comm, ierror )
  integer, intent(in) :: n
  integer, intent(in) :: data1
  integer, intent(out) :: data2
  integer, intent(in) :: op, comm, dtype
  integer, intent(out) :: ierror
  ierror = MPI_SUCCESS
  data2 = data1
  end subroutine

  subroutine mpi_allreduce_r0 (data1, data2, n, dtype, op, comm, ierror )
  use eirmod_precision
  integer, intent(in) :: n
  real(dp), intent(in) :: data1
  real(dp), intent(out) :: data2
  integer, intent(in) :: op, comm, dtype
  integer, intent(out) :: ierror
  ierror = MPI_SUCCESS
  data2 = data1
  end subroutine
#endif

#if MPI_VERSION < 3 || !defined(USE_MPI)
#ifndef GFORTRAN 
  subroutine mpi_ireduce_i0_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    use eirmod_comprt
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    request = 0
    write(iunout,*) 'Error MPI_IREDUCE called, but it is not implemented'
    call mpi_abort(MPI_COMM_WORLD, -1, ierror)
  end subroutine

  subroutine mpi_ireduce_i0_l1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    use eirmod_comprt
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    logical, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    request = 0
    write(iunout,*) 'Error MPI_IREDUCE called, but it is not implemented'
    call mpi_abort(MPI_COMM_WORLD, -1, ierror)
  end subroutine

  subroutine mpi_ireduce_r1_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    use eirmod_comprt
    implicit none
    integer, intent(in) :: n, comm
    double precision, dimension(:), intent(in) :: data1
    double precision, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    request = 0
    write(iunout,*) 'Error MPI_IREDUCE called, but it is not implemented'
    call mpi_abort(MPI_COMM_WORLD, -1, ierror)
  end subroutine

  subroutine mpi_ireduce_l1_l1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    use eirmod_comprt
    implicit none
    integer, intent(in) :: n, comm
    logical, dimension(:), intent(in) :: data1
    logical, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    request = 0
    write(iunout,*) 'Error MPI_IREDUCE called, but it is not implemented'
    call mpi_abort(MPI_COMM_WORLD, -1, ierror)
  end subroutine
#endif
#endif

#ifndef USE_MPI
  subroutine mpi_reduce_i0_i0 ( data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    integer, intent(out) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( datatype == mpi_integer ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    data2 = data1
  end subroutine

  subroutine mpi_reduce_i1_i1 ( data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1(n)
    integer, intent(out) :: data2(n)
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( datatype == mpi_integer ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    data2 = data1
  end subroutine

  subroutine mpi_reduce_L1_L1 ( data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    logical, intent(in) :: data1(n)
    logical, intent(out) :: data2(n)
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( datatype == mpi_logical ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    data2 = data1
  end subroutine

  subroutine mpi_reduce_L2_L1 ( data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    logical, intent(in) :: data1(:,:)
    logical, intent(out) :: data2(n)
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( datatype == mpi_logical ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    data2(1:n) = data1(1:n,1)
  end subroutine

#ifdef GFORTRAN
  subroutine mpi_ireduce_i0_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    if ( datatype == mpi_integer ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    request = 0
    data2 = data1
  end subroutine

  subroutine mpi_ireduce_i0_l1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    logical, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    if ( datatype == mpi_integer ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    request = 0
    if ( data1 == 0 ) then
      data2 = .false.
    else if ( data1 == 1 ) then
      data2 = .true.
    else
      ierror = MPI_FAILURE
    end if
  end subroutine
  
  subroutine mpi_ireduce_r1_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    double precision, dimension(:), intent(in) :: data1
    double precision, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    if ( datatype == mpi_double_precision ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    request = 0
    data2 = data1
  end subroutine
  
  subroutine mpi_ireduce_l1_l1 (data1, data2, n, datatype, operation, receiver, &
    comm, request, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    logical, dimension(:), intent(in) :: data1
    logical, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: request, ierror
    if ( datatype == mpi_logical ) then
      ierror = mpi_success
    else
      ierror = MPI_FAILURE
    end if
    request = 0
    data2 = data1
  end subroutine
#endif

  subroutine mpi_reduce_inplace_r0 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( data1 == MPI_IN_PLACE ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_reduce_inplace_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( data1 == MPI_IN_PLACE ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_reduce_inplace_r2 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( data1 == MPI_IN_PLACE ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
  end subroutine
  subroutine mpi_reduce_inplace_r3 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:,:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( data1 == MPI_IN_PLACE ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_reduce_inplace_r4 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    ! first argument is MPI_IN_PLACE flag
    implicit none
    integer, intent(in) :: n, comm
    integer, intent(in) :: data1
    double precision, dimension(:,:,:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( data1 == MPI_IN_PLACE ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_reduce_r0_r0 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    double precision, intent(in) :: data1 ! buffers should not overlap
    double precision, intent(out) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
    data2 = data1
  end subroutine

  subroutine mpi_reduce_r1_r1 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    double precision, intent(in), dimension(:) :: data1 ! buffers should not overlap
    double precision, intent(out), dimension(:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( sum(abs(shape(data1)-shape(data2)))==0 ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
    data2 = data1
  end subroutine

  subroutine mpi_reduce_r2_r2 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    double precision, intent(in), dimension(:,:) :: data1 ! buffers should not overlap
    double precision, intent(out), dimension(:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             (  sum(abs(shape(data1)-shape(data2)))==0 ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
    data2 = data1
  end subroutine

  subroutine mpi_reduce_r3_r3 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    double precision, intent(in), dimension(:,:,:) :: data1 ! buffers should not overlap
    double precision, intent(out), dimension(:,:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( sum(abs(shape(data1)-shape(data2)))==0 ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
    data2 = data1
  end subroutine

  subroutine mpi_reduce_r4_r4 (data1, data2, n, datatype, operation, receiver, &
    comm, ierror )
    implicit none
    integer, intent(in) :: n, comm
    double precision, intent(in), dimension(:,:,:,:) :: data1 ! buffers should not overlap
    double precision, intent(out), dimension(:,:,:,:) :: data2
    integer, intent(in) :: datatype, operation, receiver
    integer, intent(out) :: ierror
    if ( ( datatype == MPI_DOUBLE_PRECISION ) .and. &
             ( sum(abs(shape(data1)-shape(data2)))==0 ) ) then
      ierror = MPI_SUCCESS
    else
      ierror = MPI_FAILURE
    endif
    data2 = data1
  end subroutine

  subroutine mpi_gather_i0_i1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
    root, comm, ierr)
    implicit none
    integer, intent(in) :: sendbuf
    integer, intent(in) :: sendcount, sendtype, recvtype
    integer, intent(in) :: recvcount
    integer, intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    ierr = mpi_success
    if(sendbuf.ne.MPI_IN_PLACE) then
      recvbuf(1) = sendbuf
    endif
  end subroutine

  subroutine mpi_gather_i1_i1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
    root, comm, ierr)
    implicit none
    integer, intent(in) :: sendcount
    integer, intent(in) :: sendbuf(sendcount)
    integer, intent(in) :: sendtype, recvtype
    integer, intent(in) :: recvcount
    integer, intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    if (sendcount.eq.recvcount) then
      ierr = mpi_success
      recvbuf(1:recvcount) = sendbuf(1:sendcount)
    else
      ierr = mpi_failure
    endif
  end subroutine

  subroutine mpi_gather_r0_r0(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
    root, comm, ierr)
    implicit none
    integer, intent(in) :: sendcount
    double precision, intent(in) :: sendbuf
    integer, intent(in) :: sendtype, recvtype
    integer, intent(in) :: recvcount
    double precision, intent(out) :: recvbuf
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    if (sendcount.eq.recvcount) then
      ierr = mpi_success
      recvbuf = sendbuf
    else
      ierr = mpi_failure
    endif
  end subroutine

  subroutine mpi_gather_r0_r1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
    root, comm, ierr)
    implicit none
    integer, intent(in) :: sendcount
    double precision, intent(in) :: sendbuf
    integer, intent(in) :: sendtype, recvtype
    integer, intent(in) :: recvcount
    double precision, intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    ierr = mpi_success
    recvbuf(1:recvcount) = sendbuf
  end subroutine

  subroutine mpi_gather_r1_r1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
    root, comm, ierr)
    implicit none
    integer, intent(in) :: sendcount
    double precision, intent(in) :: sendbuf(sendcount)
    integer, intent(in) :: sendtype, recvtype
    integer, intent(in) :: recvcount
    double precision, intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    if (sendcount.eq.recvcount) then
      ierr = mpi_success
      recvbuf(1:recvcount) = sendbuf(1:sendcount)
    else
      ierr = mpi_failure
    endif
  end subroutine

  subroutine mpi_comm_group(comm, group, ierr)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: group
    integer, intent(out) :: ierr
    group = comm
    ierr = mpi_success
  end subroutine

  subroutine mpi_comm_split(comm, color, key, newcomm, ierr)
    implicit none
    integer, intent(in) :: comm, color, key
    integer, intent(out) :: newcomm
    integer, intent(out) :: ierr
    newcomm = comm
    ierr = mpi_success
  end subroutine

  subroutine mpi_comm_create(comm, group, newcomm, ierr)
   implicit none
   integer, intent(in) :: comm, group
   integer, intent(out) :: newcomm, ierr
   newcomm = group
   ierr = mpi_success
  end subroutine

  subroutine mpi_comm_free(comm, ierr)
   implicit none
   integer, intent(inout) :: comm
   integer, intent(out) :: ierr
   comm = mpi_comm_null
   ierr = mpi_success
  end subroutine

  subroutine mpi_comm_rank (comm, me, ierr)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: me, ierr
    me = 0
    ierr = MPI_SUCCESS
  end subroutine

  subroutine mpi_comm_size (comm, size, ierr)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: size, ierr
    if (comm.ne.MPI_COMM_NULL) then
      size = 1
    else
      size = 0
    endif
    ierr = MPI_SUCCESS
  end subroutine

  subroutine mpi_get_version(major, minor, ierr)
    implicit none
    integer, intent(out) :: major, minor, ierr
    major = MPI_VERSION
    minor = MPI_SUBVERSION
    ierr = MPI_SUCCESS
  end subroutine

  subroutine mpi_group_size(group, gsize, ierr)
    implicit none
    integer, intent(in) :: group
    integer, intent(out) :: gsize, ierr
    if (group == mpi_group_empty) then
      gsize = 0
    else
      gsize = 1
    endif
    ierr = mpi_success
  end subroutine

  subroutine mpi_group_rank(group, grank, ierr)
    implicit none
    integer, intent(in) :: group
    integer, intent(out) :: grank, ierr
    grank = 0
    ierr = mpi_success
  end subroutine

  subroutine mpi_group_incl(group, n, ranks, newgroup, ierr)
    implicit none
    integer, intent(in) :: group
    integer, intent(in) :: n
    integer, intent(in), dimension(n) :: ranks
    integer, intent(out) :: newgroup, ierr
    ierr = mpi_success
    if (n<0 .or. group == mpi_group_null) then
      newgroup = mpi_group_empty
    else
      newgroup = group
    endif
  end subroutine

  subroutine mpi_group_free(group, ierr)
    implicit none
    integer, intent(inout) :: group
    integer, intent(out) :: ierr
    group = mpi_group_null
    ierr = mpi_success
  end subroutine

  subroutine mpi_init(ierr)
    implicit none
    integer, intent(out) :: ierr
    ierr = mpi_success
  end subroutine

  subroutine mpi_barrier(comm, ierr)
    implicit none
    integer, intent(in) :: comm
    integer, intent(out) :: ierr
    ierr = mpi_success
  end subroutine

  subroutine mpi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr
    ierr = mpi_success
  end subroutine

  subroutine mpi_scatter_i0_i0(sendbuf, sendcount, sendtype,  &
             recvbuf, recvcount, recvtype, root, comm, ierr)
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount, recvcount
    integer, intent(in) :: sendbuf
    integer, intent(out) :: recvbuf
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    if (sendcount == recvcount) then
      ierr = MPI_SUCCESS
      recvbuf = sendbuf
    else
      write (iunout,*) 'MPI SCATTER: buffer size missmatch: ', &
       sendcount, recvcount
      call eirene_exit_own(1)
    endif
  end subroutine

  subroutine mpi_scatter_i1_i0(sendbuf, sendcount, sendtype,  &
             recvbuf, recvcount, recvtype, root, comm, ierr)
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount, recvcount
    integer, intent(in) :: sendbuf(sendcount)
    integer, intent(out) :: recvbuf
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    if (sendcount == recvcount) then
      ierr = MPI_SUCCESS
      recvbuf = sendbuf(1)
    else
      write (iunout,*) 'MPI SCATTER: buffer size missmatch: ', &
       sendcount, recvcount
      call eirene_exit_own(1)
    endif
  end subroutine

  subroutine mpi_scatter_i1_i1(sendbuf, sendcount, sendtype,  &
             recvbuf, recvcount, recvtype, root, comm, ierr)
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount, recvcount
    integer, intent(in) :: sendbuf(sendcount)
    integer, intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    if (sendcount == recvcount) then
      ierr = MPI_SUCCESS
      recvbuf(1:sendcount) = sendbuf(1:sendcount)
    else
      write (iunout,*) 'MPI SCATTER: buffer size missmatch: ', &
       sendcount, recvcount
      call eirene_exit_own(1)
    endif
  end subroutine

  subroutine mpi_gatherv_r2_r2 (sendbuf, sendcount, sendtype,  &
             recvbuf, recvcount, recvdistrib, recvtype, root, comm, ierr)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount, recvcount(*)
    real(dp), intent(in) :: sendbuf(:,:)
    integer, intent(in) :: recvdistrib(*)
    real(dp), intent(out) :: recvbuf(:,:)
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    integer :: i
    ierr = MPI_SUCCESS
    do i=1, sendcount
      recvbuf(recvdistrib(i),1) = sendbuf(i,1)
    end do
  end subroutine

  subroutine mpi_gatherv_i2_i2 (sendbuf, sendcount, sendtype,  &
             recvbuf, recvcount, recvdistrib, recvtype, root, comm, ierr)
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount, recvcount(*)
    integer, intent(in) :: sendbuf(:,:)
    integer, intent(in) :: recvdistrib(*)
    integer, intent(out) :: recvbuf(:,:)
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    integer :: i
    ierr = MPI_SUCCESS
    do i=1, sendcount
      recvbuf(recvdistrib(i),1) = sendbuf(i,1)
    end do
  end subroutine

  subroutine mpi_scatterv (sendbuf, sendcount, senddistrib, sendtype,  &
             recvbuf, recvcount, recvtype, root, comm, ierr)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(in) :: sendcount(*), recvcount
    real(dp), intent(in) :: sendbuf(*)
    integer, intent(in) :: senddistrib(*)
    real(dp), intent(out) :: recvbuf(recvcount)
    integer, intent(in) :: sendtype, recvtype, root, comm
    integer, intent(out) :: ierr
    integer :: i
    ierr = MPI_SUCCESS
    do i=1, recvcount
      recvbuf(i) = sendbuf(senddistrib(i))
    end do
  end subroutine

  subroutine mpi_send_dum_a0 (buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_a1 (buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:)
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_a2(buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_a3(buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:)
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_i0 (buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    integer :: buffer
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_i1 (buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:)
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_send_dum_i2 (buffer,cnt,datatype,dest,tag,comm,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,dest,tag,comm
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI send called in a serial code. No data is sent.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_a0(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    integer, intent(out) :: ier
    real(dp) :: buffer
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_a1(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:)
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_a2(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_a3(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    real(dp) :: buffer(:,:,:)
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_i0(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    integer, intent(out) :: ier
    integer :: buffer
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_i1(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:)
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_recv_dum_i2(buffer,cnt,datatype,source,tag,comm,st,ier)
    use eirmod_precision
    use eirmod_comprt
    implicit none
    integer, intent(out) :: ier
    integer :: buffer(:,:)
    integer, intent(in) :: cnt,datatype,source,tag,comm
    integer :: st
    ier=MPI_SUCCESS
    write(iunout,*) &
     'Warning: MPI recv called in a serial code. No data is received.'
    write(iunout,*) &
     '         Results are only correct if send and recv buffers are the same'
  end subroutine

  subroutine mpi_waitall_i0(n, requests, status, ierr)
  ! status must be MPI_STATUSES_IGNORE in this case
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: requests(n)
    integer, intent(in) :: status
    integer, intent(out) :: ierr
    requests = MPI_REQUEST_NULL
    if (status == MPI_STATUSES_IGNORE) then
      ierr = MPI_SUCCESS
    else
      ierr = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_waitall_i1(n, requests, statuses, ierr)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: requests(n)
    integer, dimension(MPI_STATUS_SIZE, n) :: statuses
    integer, intent(out) :: ierr
    requests = MPI_REQUEST_NULL
    ierr = MPI_SUCCESS
  end subroutine

  subroutine mpi_testall_i0(n, requests, flag, status, ierr)
  ! status must be MPI_STATUSES_IGNORE in this case
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: requests(n)
    logical, intent(out) :: flag
    integer, intent(in) :: status
    integer, intent(out) :: ierr
    requests = MPI_REQUEST_NULL
    flag = .true.
    if (status == MPI_STATUSES_IGNORE) then
      ierr = MPI_SUCCESS
    else
      ierr = MPI_FAILURE
    endif
  end subroutine

  subroutine mpi_testall_i1(n, requests, flag, statuses, ierr)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: requests(n)
    logical, intent(out) ::  flag
    integer, dimension(MPI_STATUS_SIZE, n) :: statuses
    integer, intent(out) :: ierr
    flag = .true.
    requests = MPI_REQUEST_NULL
    ierr = MPI_SUCCESS
  end subroutine
#endif
end module

!!!Local Variables:
!!! mode: f90
!!! End:
