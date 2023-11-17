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

      TREEPATH = get_solpstop()
      IF (INDEX(TREEPATH,' ').NE.1) THEN
        I4 = SCAN(TREEPATH,' ')-1
        DBFNAME(1:I4+16) = TREEPATH(1:I4)//'/modules/Eirene/'
        DBFNAME(I4+17:I4+16+IEND-IANF+1) = ZEILE(IANF:IEND)
        inquire(file=DBFNAME,exist=ex)
        if (.not.ex) then
          DBFNAME(1:I4+1) = TREEPATH(1:I4)//'/'
          DBFNAME(I4+2:I4+1+IEND-IANF+1) = ZEILE(IANF:IEND)
          DBFNAME(I4+1+IEND-IANF+2:400) =
     .         REPEAT(' ',400-(I4+1+IEND-IANF+2)+1)
        end if
      END IF

      return
      end subroutine eirene_filepath_usr
