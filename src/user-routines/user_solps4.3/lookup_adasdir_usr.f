      subroutine eirene_lookup_adasdir_usr(DSN, FOUND, 
     .                                     REAC, ELNAME, BUNDLING)
      use eirmod_cinit
      implicit none
      character(*), intent(inout) :: dsn, reac
      logical, intent(inout) :: found
      character(*), optional :: ELNAME, BUNDLING

      found = .false.

      return
      end subroutine eirene_lookup_adasdir_usr
