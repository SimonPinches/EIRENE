      subroutine eirene_lookup_adasdir_usr(DSN, FOUND, 
     .                                     REAC, ELNAME, BUNDLING)
      use eirmod_cinit
      implicit none
      character(*),intent(inout) :: DSN
      character(*), intent(in) :: REAC
      logical, intent(inout) :: found
      character(*), intent(in), optional :: ELNAME, BUNDLING

      found = .false.

      return
      end subroutine eirene_lookup_adasdir_usr
