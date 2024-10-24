
      MODULE EIRMOD_JSON

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD, ONLY: NLIM, NLIMPS,
     .                          NATM, NMOL, NION, NPHOT, NPLS, NSPZ,
     .                          NSTRA
!pb      USE EIRMOD_COMPRT, ONLY: IUNOUT

      use json_module           !IGNORE
!cym/cpg keep original names
!     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck
!cym, for reference : this would be equivalent to the initial code; works too for this file
!      use json_kinds, only: lk, rk, ik, ck

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: eirene_init_input_blocks, jval_array,
     .          eirene_alloc_json_arrays, eirene_dealloc_json_arrays,
     .          S_stack, STRING_STACK,
     .          eirene_push_string_stack, eirene_copy_addsrf,
     p          TRANSFORM, TRANSFORM_QUEUE

      integer, public, parameter :: nblks=15
      character(20), public, parameter :: blknam(0:nblks)=
     .            (/ "HEADER              ",
     1               "GENERAL_DATA        ",
     2               "STANDARD_MESH       ",
     3               "NONDEF_STD_SURFACES ",
     4               "ADDITIONAL_SURFACES ",
     5               "PHYSICS_MODEL       ",
     6               "REFLECTION_MODELS   ",
     7               "SOURCES             ",
     8               "SPECIFIC_ZONES      ",
     9               "STATISTICS          ",
     .               "ADDITIONAL_TALLIES  ",
     1               "OUTPUT              ",
     2               "DIAGNOSTICS         ",
     3               "TIMEDEPENDENT_MODE  ",
     4               "INTERFACING         ",
     5               "MPI_INFORMATION     " /)
      logical, public, parameter :: blk_required(0:nblks)=
     >     (/.true.,  .true.,  .true.,  .true.,  .true.,
     >       .true.,  .true.,  .true.,  .false., .false.,
     >       .false., .true.,  .false., .false., .true.,
     >       .false. /)

      type jval_array
        type(json_value), pointer :: p
      end type jval_array

      ! array of json trees for reading from different files
      type(json_core), public, save :: jtrees(0:nblks)
      ! pointer to first item of jtree
      type(jval_array), public, save :: ptree(0:nblks)
      ! pointer to the branch of a block in a jtree
      type(jval_array), public, save :: blks(0:nblks)
      ! number of jtrees that have been read
      integer, public, save :: njs=-1,
      ! number of the jtree containing the data of a block
     .                         itree_num(0:nblks)


!     variables stored for output to json file

! arrays for json output
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: SORCOS_IN(:), SORMAX_IN(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R P1_IN(:,:), P2_IN(:,:), P3_IN(:,:), P4_IN(:,:), P5_IN(:,:),
     R P6_IN(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R A0LM_IN(:), A1LM_IN(:), A2LM_IN(:), A3LM_IN(:), A4LM_IN(:),
     R A5LM_IN(:), A6LM_IN(:), A7LM_IN(:), A8LM_IN(:), A9LM_IN(:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R XLIMS1_IN(:,:), YLIMS1_IN(:,:), ZLIMS1_IN(:,:),
     R XLIMS2_IN(:,:), YLIMS2_IN(:,:), ZLIMS2_IN(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R ALIMS_IN(:,:), XLIMS_IN(:,:), YLIMS_IN(:,:), ZLIMS_IN(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R ALIMS0_IN(:,:), XLIMS3_IN(:,:), YLIMS3_IN(:,:), ZLIMS3_IN(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: RLB_IN(:)
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I ILCOL_IN(:), LCHSPNWL_IN(:,:)

      REAl(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R DATD_IN(:), DMLD_IN(:), DIOD_IN(:), DPLD_IN(:), DPHD_IN(:)

! storing original INDPRO data for json output
      INTEGER, PUBLIC, SAVE :: INDPRO_IN(12)

! arrays for storing numbers of tallies explicitly switched on or off
! allocated in subroutine input
      INTEGER, PUBLIC, ALLOCATABLE, SAVE :: ITLVOUT(:), ITLSOUT(:)
      INTEGER, PUBLIC, SAVE :: NTLVOUT, NTLSOUT

      INTEGER, PUBLIC, SAVE :: NR1ST_IN, NP2ND_IN, NT3RD_IN, NRTAL_IN,
     .                         NOPTIM_IN, NSMSTRA_IN, NSTRAI_IN,
     .                         NTIME_IN
      LOGICAL, PUBLIC, SAVE :: NLPOL_IN, NLPLG_IN, NLFEM_IN

      LOGICAL, PUBLIC, SAVE :: LDEF_LINES, LDEF_TIME_HORIZON


      TYPE TRANSFORM
        INTEGER :: ITINI, ITEND
        REAL(DP) :: XLCOR, YLCOR, ZLCOR,
     R              XLREF, YLREF, ZLREF,
     R              XLROT, YLROT, ZLROT, ALROT
        TYPE(TRANSFORM), POINTER :: NEXT
      END TYPE TRANSFORM

      TYPE TRANSFORM_QUEUE
        TYPE(TRANSFORM), POINTER :: HEAD, LAST
      END TYPE TRANSFORM_QUEUE

      TYPE S_STACK
#ifdef LEGACYCOMP
        CHARACTER, dimension(:), allocatable :: STRING
#else
        CHARACTER(len=:), allocatable :: STRING
#endif
        TYPE(S_STACK), POINTER :: NEXT
      END TYPE S_STACK

      TYPE STRING_STACK
        TYPE(S_STACK), POINTER :: HEAD, LAST
      END TYPE STRING_STACK

      TYPE(STRING_STACK), PUBLIC, SAVE :: CM_STACK,  CH0_STACK,
     >                                    CRS_STACK
      TYPE(STRING_STACK), PUBLIC, ALLOCATABLE, SAVE ::
     >  CH1_STACK(:), CH2_STACK(:), CH3_STACK(:)

      TYPE(TRANSFORM_QUEUE), PUBLIC, ALLOCATABLE, SAVE :: TRNSFRM(:)

      contains


!******************************************************************************

      subroutine eirene_start_json(fname,iblk,iunout)

      implicit none
      character(*), intent(in) :: fname
      integer, intent(in) :: iblk, iunout
      logical :: status_ok
      character(kind=json_CK,len=:), allocatable :: error_msg
      external :: eirene_exit_own

      njs = njs + 1

!      call jtrees(njs)%initialize(compact_reals=.true.,real_format='ES',
!     .          case_sensitive_keys=.false.,compress_vectors=.true.)
      call jtrees(njs)%initialize(compact_reals=.true.,real_format='ES',
     .          case_sensitive_keys=.false.,compress_vectors=.false.)

      call jtrees(njs)%parse(file=fname, p=ptree(njs)%p)
      call jtrees(njs)%check_for_errors(status_ok,error_msg)
!      if (jtrees(njs)%failed()) then
      if (.not. status_ok) then
         write (iunout,*) ' EIRENE INPUT FILE ',fname,
     .                    ' COULD NOT BE READ '
         write(*,*) 'Error: '//error_msg
         call eirene_exit_own(1)
      end if

      itree_num(iblk) = njs

      return
      end subroutine eirene_start_json


!******************************************************************************


      subroutine eirene_init_input_blocks (inname,iunout)

      IMPLICIT NONE
      character(*), intent(in) :: inname
      integer, intent(in) :: iunout
      integer :: i
      character(200) :: fname
!cym/cpg ck -> json_ck
      character(kind=json_CK,len=:), allocatable :: incname
!cym - avoid type mismatch
      character(kind=json_CK,len=1), parameter :: space=' '
!cym
      integer :: ind, j
      logical :: found
      external :: eirene_leer, eirene_exit_own

      call eirene_start_json(inname,0,iunout)

!  find starting points for reading the blocks of the EIRENE input file
      do i = 0, nblks

!  check for reading of a separate file for input block I
        call jtrees(0)%get(ptree(0)%p,trim(blknam(i))//'.INCLUDE',
     .                     incname,found)

        if (found) then
!  separate file found, get filename
!cym gfortran reports type mismatch
!cym          ind = index(incname,' ')
          ind = index(incname,space)
!cym
          if (ind == 0) then
            fname = incname
          else
            fname = incname(1:ind-1)
          end if
          deallocate(incname)

          write (iunout,'(A,A)')
     .       'EXTERNAL FILE FOUND FOR INPUT BLOCK: ',blknam(i)
          write (iunout,*) 'FILE = ',trim(fname)
          call eirene_leer(1)

!  initialize separate file
          call eirene_start_json(fname,i,iunout)
        end if

! index of the jtree that holds the data of the input block
        j = itree_num(i)
        call jtrees(j)%get_child(ptree(j)%p,trim(blknam(i)),
     .                           blks(i)%p,found)
        if (.not.found) then
          if (blk_required(i)) then
            write (iunout,*) ' EIRENE INPUT BLOCK ',blknam(i),
     .                       ' COULD NOT BE FOUND IN JSON FILE '
            call eirene_exit_own(1)
          else
            nullify(blks(i)%p)
          end if
        end if

      end do

      return

      end subroutine eirene_init_input_blocks


!******************************************************************************

      subroutine eirene_alloc_json_arrays

      ALLOCATE(SORCOS_IN(NSTRA))
      ALLOCATE(SORMAX_IN(NSTRA))

      ALLOCATE(RLB_IN(NLIM))

      ALLOCATE(P1_IN(3,NLIM))
      ALLOCATE(P2_IN(3,NLIM))
      ALLOCATE(P3_IN(3,NLIM))
      ALLOCATE(P4_IN(3,NLIM))
      ALLOCATE(P5_IN(3,NLIM))
      ALLOCATE(P6_IN(3,NLIM))

      ALLOCATE(A0LM_IN(NLIM))
      ALLOCATE(A1LM_IN(NLIM))
      ALLOCATE(A2LM_IN(NLIM))
      ALLOCATE(A3LM_IN(NLIM))
      ALLOCATE(A4LM_IN(NLIM))
      ALLOCATE(A5LM_IN(NLIM))
      ALLOCATE(A6LM_IN(NLIM))
      ALLOCATE(A7LM_IN(NLIM))
      ALLOCATE(A8LM_IN(NLIM))
      ALLOCATE(A9LM_IN(NLIM))

      ALLOCATE(XLIMS1_IN(9,NLIM))
      ALLOCATE(YLIMS1_IN(9,NLIM))
      ALLOCATE(ZLIMS1_IN(9,NLIM))
      ALLOCATE(XLIMS2_IN(9,NLIM))
      ALLOCATE(YLIMS2_IN(9,NLIM))
      ALLOCATE(ZLIMS2_IN(9,NLIM))
      ALLOCATE(ALIMS_IN(9,NLIM))
      ALLOCATE(XLIMS_IN(9,NLIM))
      ALLOCATE(YLIMS_IN(9,NLIM))
      ALLOCATE(ZLIMS_IN(9,NLIM))
      ALLOCATE(ALIMS0_IN(9,NLIM))
      ALLOCATE(XLIMS3_IN(9,NLIM))
      ALLOCATE(YLIMS3_IN(9,NLIM))
      ALLOCATE(ZLIMS3_IN(9,NLIM))

      ALLOCATE(ILCOL_IN(NLIMPS))
      ALLOCATE(LCHSPNWL_IN(NSPZ,0:NLIMPS))

      ALLOCATE(DATD_IN(NATM))
      ALLOCATE(DMLD_IN(NMOL))
      ALLOCATE(DIOD_IN(NION))
      ALLOCATE(DPLD_IN(NPLS))
      ALLOCATE(DPHD_IN(NPHOT))

      ALLOCATE (CH1_STACK(NLIM))
      ALLOCATE (CH2_STACK(NLIM))

      ALLOCATE (TRNSFRM(NLIM))

      call eirene_init_json_arrays

      end subroutine eirene_alloc_json_arrays

!******************************************************************************

      subroutine eirene_dealloc_json_arrays
      INTEGER :: I

      DEALLOCATE(SORCOS_IN)
      DEALLOCATE(SORMAX_IN)

      DEALLOCATE(RLB_IN)

      DEALLOCATE(P1_IN)
      DEALLOCATE(P2_IN)
      DEALLOCATE(P3_IN)
      DEALLOCATE(P4_IN)
      DEALLOCATE(P5_IN)
      DEALLOCATE(P6_IN)

      DEALLOCATE(A0LM_IN)
      DEALLOCATE(A1LM_IN)
      DEALLOCATE(A2LM_IN)
      DEALLOCATE(A3LM_IN)
      DEALLOCATE(A4LM_IN)
      DEALLOCATE(A5LM_IN)
      DEALLOCATE(A6LM_IN)
      DEALLOCATE(A7LM_IN)
      DEALLOCATE(A8LM_IN)
      DEALLOCATE(A9LM_IN)

      DEALLOCATE(XLIMS1_IN)
      DEALLOCATE(YLIMS1_IN)
      DEALLOCATE(ZLIMS1_IN)
      DEALLOCATE(XLIMS2_IN)
      DEALLOCATE(YLIMS2_IN)
      DEALLOCATE(ZLIMS2_IN)
      DEALLOCATE(ALIMS_IN)
      DEALLOCATE(XLIMS_IN)
      DEALLOCATE(YLIMS_IN)
      DEALLOCATE(ZLIMS_IN)
      DEALLOCATE(ALIMS0_IN)
      DEALLOCATE(XLIMS3_IN)
      DEALLOCATE(YLIMS3_IN)
      DEALLOCATE(ZLIMS3_IN)

      DEALLOCATE(ILCOL_IN)
      DEALLOCATE(LCHSPNWL_IN)

      DEALLOCATE(DATD_IN)
      DEALLOCATE(DMLD_IN)
      DEALLOCATE(DIOD_IN)
      DEALLOCATE(DPLD_IN)
      DEALLOCATE(DPHD_IN)

      call eirene_resolve_string_stack(cm_stack)
      call eirene_resolve_string_stack(ch0_stack)
      call eirene_resolve_string_stack(crs_stack)

      DO I=1,NLIM
        call eirene_resolve_string_stack(ch1_stack(i))
        call eirene_resolve_string_stack(ch2_stack(i))
      END DO

      DEALLOCATE (CH1_STACK)
      DEALLOCATE (CH2_STACK)

      DEALLOCATE (TRNSFRM)

      end subroutine eirene_dealloc_json_arrays

!******************************************************************************

      subroutine eirene_resolve_string_stack (stack)

      type(string_stack), intent(inout) :: stack
      TYPE(S_STACK), POINTER :: CUR

      CUR => STACK%HEAD
      DO WHILE (ASSOCIATED(CUR))
        STACK%HEAD => CUR%NEXT
        DEALLOCATE (CUR)
        CUR => STACK%HEAD
      END DO
      nullify(stack%head)
      nullify(stack%last)

      end subroutine eirene_resolve_string_stack

!******************************************************************************

      subroutine eirene_init_json_arrays
      integer :: i

      SORCOS_IN = 0._DP
      SORMAX_IN = 0._DP

      RLB_IN = 0._DP

      P1_IN = 0._DP
      P2_IN = 0._DP
      P3_IN = 0._DP
      P4_IN = 0._DP
      P5_IN = 0._DP
      P6_IN = 0._DP

      A0LM_IN = 0._DP
      A1LM_IN = 0._DP
      A2LM_IN = 0._DP
      A3LM_IN = 0._DP
      A4LM_IN = 0._DP
      A5LM_IN = 0._DP
      A6LM_IN = 0._DP
      A7LM_IN = 0._DP
      A8LM_IN = 0._DP
      A9LM_IN = 0._DP

      XLIMS1_IN = 0._DP
      YLIMS1_IN = 0._DP
      ZLIMS1_IN = 0._DP
      XLIMS2_IN = 0._DP
      YLIMS2_IN = 0._DP
      ZLIMS2_IN = 0._DP
      ALIMS_IN = 0._DP
      XLIMS_IN = 0._DP
      YLIMS_IN = 0._DP
      ZLIMS_IN = 0._DP
      ALIMS0_IN = 0._DP
      XLIMS3_IN = 0._DP
      YLIMS3_IN = 0._DP
      ZLIMS3_IN = 0._DP

      ILCOL_IN = 0
      LCHSPNWL_IN = 0

      DATD_IN = 0._DP
      DMLD_IN = 0._DP
      DIOD_IN = 0._DP
      DPLD_IN = 0._DP
      DPHD_IN = 0._DP

      NULLIFY (CM_STACK%HEAD)

      NULLIFY(CH0_STACK%HEAD)

      NULLIFY(CRS_STACK%HEAD)

      DO I=1,NLIM
        NULLIFY(CH1_STACK(I)%HEAD)
        NULLIFY(CH1_STACK(I)%LAST)
        NULLIFY(CH2_STACK(I)%HEAD)
        NULLIFY(CH2_STACK(I)%LAST)
        NULLIFY(TRNSFRM(I)%HEAD)
        NULLIFY(TRNSFRM(I)%LAST)
      END DO

      end subroutine eirene_init_json_arrays

!******************************************************************************

      subroutine eirene_push_string_stack (stack, str)
!cym/cpg
!      character(*), intent(in) :: str
      character(kind=json_CK,len=*), intent(in) :: str
!cym/cpg end
      type(string_stack), intent(inout) :: stack
      type(s_stack), pointer :: new_elem

      allocate (new_elem)
      new_elem%string = str
      nullify(new_elem%next)

      if (associated(stack%head)) then
        stack%last%next => new_elem
      else
        stack%head => new_elem
      end if

      stack%last => new_elem

      END subroutine eirene_push_string_stack
!******************************************************************************

      subroutine eirene_copy_addsrf(i)
      use eirmod_cadgeo
      use eirmod_clgin
      integer, intent(in) :: i

      RLB_IN(I) = RLB(I)

      P1_IN(:,I) = P1(:,I)
      P2_IN(:,I) = P2(:,I)
      P3_IN(:,I) = P3(:,I)
      P4_IN(:,I) = P4(:,I)
      P5_IN(:,I) = P5(:,I)
      P6_IN(:,I) = P6(:,I)

      A0LM_IN(I) = A0LM(I)
      A1LM_IN(I) = A1LM(I)
      A2LM_IN(I) = A2LM(I)
      A3LM_IN(I) = A3LM(I)
      A4LM_IN(I) = A4LM(I)
      A5LM_IN(I) = A5LM(I)
      A6LM_IN(I) = A6LM(I)
      A7LM_IN(I) = A7LM(I)
      A8LM_IN(I) = A8LM(I)
      A9LM_IN(I) = A9LM(I)

      XLIMS1_IN(:,I) = XLIMS1(:,I)
      YLIMS1_IN(:,I) = YLIMS1(:,I)
      ZLIMS1_IN(:,I) = ZLIMS1(:,I)
      XLIMS2_IN(:,I) = XLIMS2(:,I)
      YLIMS2_IN(:,I) = YLIMS2(:,I)
      ZLIMS2_IN(:,I) = ZLIMS2(:,I)
      ALIMS_IN(:,I)  = ALIMS(:,I)
      XLIMS_IN(:,I)  = XLIMS(:,I)
      YLIMS_IN(:,I)  = YLIMS(:,I)
      ZLIMS_IN(:,I)  = ZLIMS(:,I)
      ALIMS0_IN(:,I) = ALIMS0(:,I)
      XLIMS3_IN(:,I) = XLIMS3(:,I)
      YLIMS3_IN(:,I) = YLIMS3(:,I)
      ZLIMS3_IN(:,I) = ZLIMS3(:,I)

      ILCOL_IN(I) = ILCOL(I)
      LCHSPNWL_IN(:,I) = LCHSPNWL(:,I)

      end subroutine eirene_copy_addsrf
!******************************************************************************

      END MODULE EIRMOD_JSON
