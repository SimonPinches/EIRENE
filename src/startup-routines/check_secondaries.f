      subroutine eirene_check_secondaries (ifrst, iscnd, ithrd, ifrth, 
     .                                     ierr)

      implicit none

      integer, intent(in) :: ifrst, iscnd, ithrd, ifrth
      integer, intent(out) :: ierr
      integer :: icount, eirene_idez
      integer :: isp(4)

      icount = 0
      ierr = 0
      
!  count number of bulk secondaries

      if (eirene_idez(ifrst,1,3) == 4) then
         icount = icount + 1
         isp(icount) = eirene_idez(ifrst,3,3)
      end if
      
      if (eirene_idez(iscnd,1,3) == 4) then
         icount = icount + 1
         isp(icount) = eirene_idez(iscnd,3,3)
      end if
      
      if (eirene_idez(ithrd,1,3) == 4) then
         icount = icount + 1
         isp(icount) = eirene_idez(ithrd,3,3)
      end if
      
      if (eirene_idez(ifrth,1,3) == 4) then
         icount = icount + 1
         isp(icount) = eirene_idez(ifrth,3,3)
      end if
      
      if (icount <= 1) return

!  more than one bulk secondaries

      if (count(isp(1:icount) == isp(1)) /= icount) ierr = 1

      return
      
      end subroutine eirene_check_secondaries
