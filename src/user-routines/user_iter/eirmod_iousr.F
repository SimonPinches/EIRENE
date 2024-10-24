      module eirmod_iousr

C  User specific read routines

      USE EIRMOD_PRECISION
      use eirmod_parmmod
      use eirmod_extrab25
      use eirmod_solps
      IMPLICIT NONE

      private

      public :: eirene_read_block_11_usr, eirene_write_block_11_usr

      interface eirene_read_block_11_usr
        procedure :: eirene_read_block_11_usr_fixed
        procedure :: eirene_read_block_11_usr_json
      end interface eirene_read_block_11_usr

      interface eirene_write_block_11_usr
        procedure :: eirene_read_block_11_usr_json
      end interface eirene_write_block_11_usr

      contains

      subroutine eirene_read_block_11_usr_fixed

      use eirmod_comprt
     , , only: iunin, iunout

      implicit none
      integer :: i, j, idumm(9)
      external :: eirene_exit_own

      backspace iunin
csw 22jul2011
      read (iunin,6666) idumm, nsrfcls
      if(abs(nsrfcls) > 0) call eirene_extrab25_srfprvsl_init()
      do j=1,abs(nsrfcls)
        read (iunin,6666) msrfcls(j)
        if(msrfcls(j).gt.0) then !{
          if(msrfcls(j).gt.msrfclsx) then !{
            write(iunout,*)
     ,          'Too many selection criteria for surface ',
     ,          'properties visualisation. j,msrfcls(j)=',j,msrfcls(j)
            CALL EIRENE_EXIT_OWN(1)
          end if !}
          read (iunin,6666) (lsrfcls(i,j),i=1,msrfcls(j))
        end if
      end do
csw
 6666 format(12i6)

      return
      end subroutine eirene_read_block_11_usr_fixed


      subroutine eirene_read_block_11_usr_json(json, me)

      use eirmod_comprt
     , , only: iunout
      use eirmod_json
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      class(json_core),intent(inout) :: json
      type(json_value),pointer :: me
      type(json_value),pointer :: pcrits, pcrit
      integer :: nsf, j
      integer, allocatable :: idum(:)
      logical :: found, foundc
      external :: eirene_exit_own

      call json%get(me,'NSRFCLS',nsrfcls,found)

      if(abs(nsrfcls) > 0) then

        call eirene_extrab25_srfprvsl_init()

        call json%get_child(me,'SURFACE_CRITERIA',pcrits,foundc)
        if (foundc) then
          call json%info(pcrits,n_children=nsf)
          if (nsf /= abs(nsrfcls)) then
            write (iunout,*)
     .      ' NUMBER OF CHILDREN IN SURFACE_CRITERIA BLOCK DOES',
     .      ' NOT MATCH NSRFCLS '
            write (iunout,*) 'NSRFCLS = ',nsrfcls
            write (iunout,*) 'NSF =   ',nsf
            call eirene_exit_own(1)
          end if

          do j=1,abs(nsrfcls)
            call json%get_child(pcrits,j,pcrit,found)
            call json%get(pcrit,'MSRFCLS', msrfcls(j),found)
            if(msrfcls(j).gt.0) then !{
              if(msrfcls(j).gt.msrfclsx) then !{
                write(iunout,*)
     ,            'Too many selection criteria for surface ',
     ,            'properties visualisation. j,msrfcls(j)=',j,msrfcls(j)
                CALL EIRENE_EXIT_OWN(1)
              end if !}
              call json%get(pcrit,'LSRFCLS',idum,found)
              lsrfcls(1:msrfcls(j),j) = idum(1:msrfcls(j))
              deallocate(idum)
            end if
          end do
        end if ! foundc
      end if

      return
      end subroutine eirene_read_block_11_usr_json


      subroutine eirene_write_block_11_usr_json(json, me)

      use eirmod_json
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      class(json_core),intent(inout) :: json
      type(json_value),pointer :: me
      type(json_value),pointer :: pcrits, pcrit
      integer :: j

      call json%add(me,'NSRFCLS',nsrfcls)

      if (abs(nsrfcls) > 0) then
        call json%create_array(pcrits,'SURFACE_CRITERIA')
        do j=1,abs(nsrfcls)
          call json%create_object(pcrit,'')
          call json%add(pcrit,'MSRFCLS',msrfcls(j))
          call json%add(pcrit,'LSRFCLS',lsrfcls(1:msrfcls(j),j))
          call json%add(pcrits,pcrit)
        end do
        call json%add(me,pcrits)
      end if

      return
      end subroutine eirene_write_block_11_usr_json

      end module eirmod_iousr


