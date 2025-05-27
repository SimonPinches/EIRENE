      subroutine eirene_filepath_usr(zeile, dbfname, ianf, iend)
      implicit none
      character(*), intent(in) :: zeile
      character(*), intent(inout) :: dbfname
      integer, intent(in) :: ianf, iend
      integer :: i4
      logical :: ex
      character(400) :: treepath
      character*256 :: get_solpstop
      external get_solpstop

      ! Check if the path is already valid
      inquire(file=DBFNAME,exist=ex)
      if (ex) return
      ! Path is not valid, but it is absolute,
      ! so no sense in prepending anything
      if (DBFNAME(1:1).eq.'/') return

      TREEPATH = get_solpstop()
      if (index(TREEPATH,' ').eq.1) return

      ! Attempt 1: Prepend $SOLPSTOP/modules/Eirene/
      I4 = scan(TREEPATH,' ')-1
      DBFNAME(1:I4+16) = TREEPATH(1:I4) // '/modules/Eirene/'
      DBFNAME(I4+17:I4+16+IEND-IANF+1) = ZEILE(IANF:IEND)
      inquire(file=DBFNAME,exist=ex)
      if (ex) return

      ! Attempt 2: Prepend $SOLPSTOP/
      DBFNAME = repeat(' ',len(DBFNAME))
      DBFNAME(1:I4+1) = TREEPATH(1:I4) // '/'
      DBFNAME(I4+2:I4+1+IEND-IANF+1) = ZEILE(IANF:IEND)

      return
      end subroutine eirene_filepath_usr
