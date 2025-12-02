      module eirmod_iousr

C  User specific read routines

      USE EIRMOD_PRECISION
      use eirmod_parmmod
      IMPLICIT NONE

      private

#ifndef NO_JSON
      public :: eirene_read_block_11_usr, eirene_write_block_11_usr

      interface eirene_read_block_11_usr
        procedure :: eirene_read_block_11_usr_fixed
        procedure :: eirene_read_block_11_usr_json
      end interface eirene_read_block_11_usr

      interface eirene_write_block_11_usr
        procedure :: eirene_write_block_11_usr_json
      end interface eirene_write_block_11_usr
#else
      public :: eirene_read_block_11_usr_fixed,
     .          eirene_write_block_11_usr
#endif

      contains

      subroutine eirene_read_block_11_usr_fixed

      use eirmod_comprt
     , , only: iunin, iunout

      implicit none
      integer :: j, idumm(9)

      return
      end subroutine eirene_read_block_11_usr_fixed

#ifndef NO_JSON
      subroutine eirene_read_block_11_usr_json(json, me)

      use eirmod_json
      use json_module           !IGNORE
!pb  .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      class(json_core) :: json
      type(json_value),pointer :: me
      type(json_value),pointer :: pcrits, pcrit
      integer :: nsf
      integer, allocatable :: idum(:)
      logical :: found, foundc

      return
      end subroutine eirene_read_block_11_usr_json


      subroutine eirene_write_block_11_usr_json(json, me)

      use eirmod_json
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      class(json_core) :: json
      type(json_value),pointer :: me
      type(json_value),pointer :: pcrits, pcrit

      return
      end subroutine eirene_write_block_11_usr_json
#else
      subroutine eirene_write_block_11_usr
      return
      end subroutine eirene_write_block_11_usr
#endif

      end module eirmod_iousr
